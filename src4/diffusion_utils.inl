/******************************************************************************
 *
 * Copyright (C) 2019-2020 by
 * The Salk Institute for Biological Studies and
 * Pittsburgh Supercomputing Center, Carnegie Mellon University
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
 * USA.
 *
******************************************************************************/

/**
 * This file is directly included into diffuse_react_event.cpp.
 * The reason why this is not a standard .cpp + .h file is to gove the compiler
 * the opportunity to inline these functions into methods of diffuse&react event.
 */
#include <vector>

#include "bng/bng.h"

#include "logging.h"

#include "diffuse_react_event.h"
#include "defines.h"
#include "world.h"
#include "partition.h"
#include "geometry.h"
#include "simulation_config.h"
#include "debug_config.h"

#include "grid_utils.inl"
#include "rxn_utils.inl"
#include "wall_utils.inl"

using namespace std;

// this makes performance worse in benchmarks, also making it configurable through API
// blocks several optimization opportunities and make causes up to 4% performance drop
//#define ENABLE_SAFE_DIFFUSION_STEP

namespace MCell {

namespace DiffusionUtils {

const double MULTISTEP_WORTHWHILE = 2.0; // TODO: make this tweakable through API, 2.0 is the value used in MCell3
const double MULTISTEP_PERCENTILE = 0.99;
const double MULTISTEP_FRACTION = 0.9;


/*************************************************************************
pick_2D_displacement:
  In: v: vector2 to store the new displacement
      scale: scale factor to apply to the displacement
      rng:
  Out: No return value.  vector is set to a random orientation and a
         distance chosen from the probability distribution of a diffusing
         2D molecule, scaled by the scaling factor.
*************************************************************************/
static void pick_surf_displacement(Vec2& v, const double scale, rng_state& rng) {
  const pos_t one_over_2_to_16th = 1.52587890625e-5f;
  Vec2 a;

  /*
   * NOTE: The below algorithm is the polar method due to Marsaglia
   * combined with a rejection method for picking uniform random
   * variates in C2.
   * Both methods are nicely described in Chapters V.4.3 and V.4.4
   * of "Non-Uniform Random Variate Generation" by Luc Devroye
   * (http://luc.devroye.org/rnbookindex.html).
   */
  pos_t f;
  do {
    unsigned int n = rng_uint(&rng);

    a.u = 2 * one_over_2_to_16th * (n & 0xFFFF) - 1;
    a.v = 2 * one_over_2_to_16th * (n >> 16) - 1;
    f = len2_squared(a);
  } while ((f < POS_EPS) || (f > 1));

  /*
   * NOTE: The scaling factor to go from a uniform to
   * a normal distribution is sqrt(-2log(f)/f).
   * However, since we use two normally distributed
   * variates to generate a normally distributed
   * 2d vector (with variance 1) we have to normalize
   * and divide by an additional factor of sqrt(2)
   * resulting in normalFactor.
   */
  pos_t normal_factor = sqrt_p(-log_p(f) / f);
  v = a * Vec2(normal_factor * (pos_t)scale);
}


static void compute_surf_displacement(
    const BNG::Species& sp,
    const double scale,
    rng_state& rng,
    Vec2& v
) {
  if (sp.can_diffuse()) {
    pick_surf_displacement(v, scale, rng);
  }
  else {
    v = Vec2(0);
  }
}


// ---------------------------------- volume mol diffusion ----------------------------------

// get displacement based on scale (related to diffusion constant) and gauss random number
static inline void pick_vol_displacement(const BNG::Species& sp, const double scale, rng_state& rng, Vec3& displacement) {

  assert(sp.can_diffuse());
  displacement.x = scale * rng_gauss(&rng) * 0.70710678118654752440;
  displacement.y = scale * rng_gauss(&rng) * 0.70710678118654752440;
  displacement.z = scale * rng_gauss(&rng) * 0.70710678118654752440;
}


#ifdef ENABLE_SAFE_DIFFUSION_STEP
/****************************************************************************
safe_diffusion_step:
  In: vm: molecule that is moving
      shead: linked list of potential collisions with molecules from the
      radial_subdivisions:
      r_step:
      x_fineparts:
      y_fineparts:
      z_fineparts:
  Out: The estimated number of diffusion steps this molecule can take before
       something interesting might happen to it, or 1.0 if something might
       happen within one timestep.
  Note: Each molecule uses its own timestep.  Only molecules that the moving
        molecule can react with directly are counted (secondary reaction
        products are ignored, so might be skipped).  "Might happen" is to
        the 99% confidence level (i.e. the distance you'd have to go before
        1% of the molecules will have gotten far enough to have a chance of
        reacting, although those 1% will probably not go in the right
        direction).  This doesn't take into account the diffusion of other
        target molecules, so it may introduce errors for clouds of molecules
        diffusing into each other from a distance.
        *FIXME*: Add a flag to make this be very conservative or to turn
        this off entirely, aside from the TIME_STEP_MAX= directive.
****************************************************************************/
// TODO: this needs to be optimized because we are losing ~7% in benchmarks
static inline double get_safe_diffusion_step(
    const Partition& p,
    const BNG::Species& sp,
    const Molecule& vm
) {
  assert(vm.is_vol());

  pos_t d2_nearmax;
  pos_t d2min = POS_FLT_GIGANTIC;
  pos_t steps;

  d2_nearmax = sp.space_step * p.config.radial_3d_step[(int)(p.config.num_radial_subdivisions * MULTISTEP_PERCENTILE)];
  d2_nearmax *= d2_nearmax;


  double subpart_side_length = p.config.subpartition_edge_length;
  IVec3 subpart_indices;
  p.get_subpart_3d_indices_from_index(vm.v.subpart_index, subpart_indices);
  Vec3 llf_subpart_boundary = p.config.partition0_llf + Vec3(subpart_indices) * Vec3(subpart_side_length);
  const Vec3& pos = vm.v.pos;

  pos_t d2;
  d2 = (pos.x - llf_subpart_boundary.x);
  d2 *= d2;
  if (d2 < d2min) {
    d2min = d2;
  }

  d2 = (llf_subpart_boundary.x + subpart_side_length - pos.x);
  d2 *= d2;
  if (d2 < d2min) {
    d2min = d2;
  }

  d2 = (pos.y - llf_subpart_boundary.y);
  d2 *= d2;
  if (d2 < d2min) {
    d2min = d2;
  }

  d2 = (llf_subpart_boundary.y + subpart_side_length - pos.y);
  d2 *= d2;
  if (d2 < d2min) {
    d2min = d2;
  }

  d2 = (pos.z - llf_subpart_boundary.z);
  d2 *= d2;
  if (d2 < d2min) {
    d2min = d2;
  }

  d2 = (llf_subpart_boundary.z + subpart_side_length - pos.z);
  d2 *= d2;
  if (d2 < d2min) {
    d2min = d2;
  }

  if (d2min < d2_nearmax) {
    return 1.0;
  }

  if (sp.has_flag(BNG::SPECIES_FLAG_CAN_VOLVOL) && !sp.has_flag(BNG::SPECIES_MOL_FLAG_CANT_INITIATE) ) {

    // get the closest distance to any molecule that can react
    const MoleculeIdsSet& sp_reactants = p.get_volume_molecule_reactants(vm.v.subpart_index, vm.species_id);
    for (molecule_id_t subpart_vm_id: sp_reactants) {
      const Molecule& subpart_vm = p.get_m(subpart_vm_id);
      assert(subpart_vm.is_vol());

      pos_t md2 = distance3_squared(vm.v.pos, subpart_vm.v.pos);

      if (md2 < d2min) {
        d2min = md2;
      }
    }
  }

  // get the closest distance to a wall in this subpart
  const WallsInSubpart& wall_ids = p.get_subpart_wall_indices(vm.v.subpart_index);
  for (wall_index_t wi: wall_ids) {
    const Wall& w = p.get_wall(wi);

    double wd2 = dot(w.normal, vm.v.pos) - w.distance_to_origin;
    wd2 *= wd2;
    if (wd2 < d2min)
      d2min = wd2;
  }


#ifdef MCELL3_4_SAFE_DIFF_STEP_RETURNS_CONSTANT
  return 1.5;
#else
  if (d2min < d2_nearmax) {
    steps = 1.0;
  }
  else {
    double steps_sq = d2min / d2_nearmax;
    if (steps_sq < MULTISTEP_WORTHWHILE * MULTISTEP_WORTHWHILE)
      steps = 1.0;
    else
      steps = sqrt_f(steps_sq);
  }
  return steps;
#endif
}
#endif


/*************************************************************************
erfcinv:

  Fast rational function approximation to inverse of the error function,
  based upon algorithm for inverse of Normal cumulative distribution
  function at http://home.online.no/~pjacklam/notes/invnorm/index.html
  by Peter J. Acklam. Accurate to about 4e-9 in absolute value.

  In: a value between 0 and 1 (not including endpoints)
  Out: the value y such that erfc(y) = input value
*************************************************************************/
static double erfcinv_f(double x) {
  /* Misc constants */
  const double tail_cutoff = 0.0485;
  const double neg_twice_log_half = 1.386294361119891;
  // const double sqrt_half_pi = 1.253314137315501;  /* For refinement */
  const double scaling_const = -0.7071067811865475;

  /* Tail numerator */
  const double tn0 = 2.938163982698783;
  const double tn1 = 4.374664141464968;
  const double tn2 = -2.549732539343734;
  const double tn3 = -2.400758277161838;
  const double tn4 = -3.223964580411365e-1;
  const double tn5 = -7.784894002430293e-3;
  /* Tail denominator */
  const double td1 = 3.754408661907416;
  const double td2 = 2.445134137142996;
  const double td3 = 3.224671290700398e-1;
  const double td4 = 7.784695709041462e-3;

  /* Central numerator */
  const double cn0 = 2.506628277459239;
  const double cn1 = -3.066479806614716e1;
  const double cn2 = 1.383577518672690e2;
  const double cn3 = -2.759285104469687e2;
  const double cn4 = 2.209460984245205e2;
  const double cn5 = -3.969683028665376e1;
  /* Central denominator */
  const double cd1 = -1.328068155288572e1;
  const double cd2 = 6.680131188771972e1;
  const double cd3 = -1.556989798598866e2;
  const double cd4 = 1.615858368580409e2;
  const double cd5 = -5.447609879822406e1;

  double p, q, r;

  if (x < tail_cutoff) {
    p = sqrt_f(-2 * log_f(x) + neg_twice_log_half);
    r = (tn0 + p * (tn1 + p * (tn2 + p * (tn3 + p * (tn4 + p * tn5))))) /
        (1.0 + p * (td1 + p * (td2 + p * (td3 + p * td4))));
  } else {
    p = 0.5 * x - 0.5;
    q = p * p;
    r = p * (cn0 + q * (cn1 + q * (cn2 + q * (cn3 + q * (cn4 + q * cn5))))) /
        (1.0 + q * (cd1 + q * (cd2 + q * (cd3 + q * (cd4 + q * cd5)))));
  }
  return scaling_const * r;
  /*
  Use the code below to refine to macine precision.  Rather slow, though.
  p = (erfc(scaling_const*r)-x)*sqrt_half_pi*exp(0.5*r*r);
  return scaling_const*(r - p/(1 + r*p/2));
  */
}


/*************************************************************************
pick_clamped_displacement:
  In: v: vector3 to store the new displacement
      vm: molecule that just came through the surface
      rng:
      radial_subdivisions:
  Out: No return value.  vector is set to a random orientation and a
         distance chosen from the probability distribution of a diffusing
         3D molecule that has come through a surface from a uniform
         concentration on the other side.
  Note: vm->previous_wall points to the wall we're coming from, and
        vm->index is the orientation we came off with
*************************************************************************/
static void pick_clamped_displacement(
    const Partition& p,
    const BNG::Species& sp,
    const Molecule& vm,
    rng_state& rng,
    Vec3& displacement) {

  const double one_over_2_to_20th = 9.5367431640625e-7;

  assert(vm.v.previous_wall_index != WALL_INDEX_INVALID);
  const Wall& w = p.get_wall(vm.v.previous_wall_index);

  uint n = rng_uint(&rng);

  /* Correct distribution along normal from surface (from lookup table) */
  double r_n = p.config.radial_2d_step[n & (p.config.num_radial_subdivisions - 1)];

  double pval = one_over_2_to_20th * ((n >> 12) + 0.5);
  double t = r_n / erfcinv_f(pval * erfc(r_n));
  Vec2 r_uv;
  pick_surf_displacement(r_uv, sqrt_f(t) * sp.space_step, rng);

  r_n *= vm.get_clamp_orientation() * sp.space_step;

  displacement = Vec3(r_n) * w.normal + Vec3(r_uv.u) * w.unit_u + Vec3(r_uv.v) * w.unit_v;
}


// - determine how far will our diffused molecule move
// - called compute_displacement in MCell3
static void compute_vol_displacement(
    const Partition& p,
    const BNG::Species& sp,
    Molecule& vm,
    double& max_time, // gets updated to t_steps
    rng_state& rng,
    Vec3& displacement,
    double& rate_factor,
    double& r_rate_factor,
    double& steps, // number of steps
    double& t_steps // time of steps
) {
  assert(max_time != 0);
  assert(vm.is_vol());

  if (vm.has_flag(MOLECULE_FLAG_ACT_CLAMPED)) {
    pick_clamped_displacement(p, sp, vm, rng, displacement);
    vm.v.previous_wall_index = WALL_INDEX_INVALID;
    // not sure why this is set, maybe because it won't be used anymore?
    vm.set_clamp_orientation(ORIENTATION_DOWN);
    vm.clear_flag(MOLECULE_FLAG_ACT_CLAMPED);

    rate_factor = 1.0;
    r_rate_factor = 1.0;
    steps = 1.0;
    t_steps = sp.time_step;
  }
  else {

  #ifdef ENABLE_SAFE_DIFFUSION_STEP
    if (max_time > MULTISTEP_WORTHWHILE) {
      // does not work correctly yet, clamping by max_time is done differently in MCell3
      steps = get_safe_diffusion_step(p, sp, vm);
    }
    else {
      steps = 1.0;
    }
  #else
    steps = 1.0;
  #endif

    t_steps = steps * sp.get_time_step();
    // clamp to max_time
    if (t_steps > max_time) {
      t_steps = max_time;
      steps = max_time / sp.get_time_step();
    }
    if (steps < EPS) {
      steps = EPS;
      t_steps = EPS * sp.get_time_step();
    }

    if (steps == 1.0) {
      pick_vol_displacement(sp, sp.get_space_step(), rng, displacement);
      r_rate_factor = rate_factor = 1.0;
    } else {
      rate_factor = sqrt_f(steps);
      r_rate_factor = 1.0 / rate_factor;
      pick_vol_displacement(sp, rate_factor * sp.get_space_step(), rng, displacement);
    }
  }

  // update max_time
  max_time = t_steps;

  p.stats.inc_diffusion_cummtime(steps);
}

// ---------------------------------- surface mol diffusion ----------------------------------



/*************************************************************************
move_sm_on_same_triangle:

  This is a helper function for diffuse_2D.

  In: world: simulation state
      sm: molecule that is moving
      new_loc: this is the location we are moving to.
      previous_box: this is the periodic box we were in previously.
      new_wall: this is the new wall we ended up on
      hd_info:
  Out: Returns true if the new location is ok and the molecule was placed
       there, otherwise returns false and the molecule must be placed
       elsewhere.
*************************************************************************/
static bool move_sm_on_same_triangle(
    Partition& p,
    Molecule& sm,
    Vec2& new_loc
) {
  Wall& wall = p.get_wall(sm.s.wall_index);
  Grid& grid = wall.grid;
  assert(grid.get_molecule_on_tile(sm.s.grid_tile_index) == sm.id);

  unsigned int new_tile_index = GridUtils::uv2grid_tile_index(new_loc, wall);

  if (new_tile_index >= grid.num_tiles) {
    mcell_internal_error("After ray_trace_2D, selected u, v coordinates "
                         "map to an out-of-bounds grid cell.  uv=(%.2f, "
                         "%.2f) sm=%d/%d",
                         new_loc.u, new_loc.v, new_tile_index, grid.num_tiles);
  }

  // We're on a new part of the grid
  molecule_id_t molecule_id_on_tile = grid.get_molecule_on_tile(new_tile_index);
  if (new_tile_index != sm.s.grid_tile_index) {
    if (molecule_id_on_tile != MOLECULE_ID_INVALID) {
      return false; /* Pick again--full here */
    }

    grid.reset_molecule_tile(sm.s.grid_tile_index);
    grid.set_molecule_tile(new_tile_index, sm.id);
    sm.s.grid_tile_index = new_tile_index;
  }

  sm.s.pos = new_loc;
  return true;
}


/*************************************************************************
move_sm_to_new_triangle:

  This is a helper function for diffuse_2D.

  In: world: simulation state
      sm: molecule that is moving
      new_loc: this is the location we are moving to.
      previous_box: this is the periodic box we were in previously.
      new_wall: this is the new wall we ended up on
      hd_info:
  Out: Returns true if the new location is ok and the molecule was placed
       there, otherwise returns false and the molecule must be placed
       elsewhere.
*************************************************************************/
static bool move_sm_to_new_triangle(
    Partition& p,
    Molecule& sm,
    Vec2& new_loc,
    const wall_index_t new_wall_index
) {
  Wall& wall = p.get_wall(sm.s.wall_index);
  Grid& grid = wall.grid;
  assert(grid.get_molecule_on_tile(sm.s.grid_tile_index) == sm.id
      && "Mapping grid tile->molecule and molecule->grid tile does not match");

  p.stats.inc_mol_moves_between_walls();

  Wall& new_wall = p.get_wall(new_wall_index);

  // No SM has been here before, so we need to make a grid on this wall.
  if (!new_wall.has_initialized_grid()) {
    new_wall.initialize_grid(p);
  }

  Grid& new_grid = new_wall.grid;

  /* Move to new tile */
  unsigned int new_tile_index = GridUtils::uv2grid_tile_index(new_loc, new_wall);

  if (new_tile_index >= new_grid.num_tiles) {
    mcell_internal_error(
        "After ray_trace_2D to a new wall, selected u, v coordinates map "
        "to an out-of-bounds grid cell.  uv=(%.2f, %.2f) sm=%d/%d",
        new_loc.u, new_loc.v, new_tile_index, new_grid.num_tiles);
  }

  molecule_id_t molecule_id_on_tile = new_grid.get_molecule_on_tile(new_tile_index);
  if (molecule_id_on_tile != MOLECULE_ID_INVALID) {
    return false; /* Pick again--full here */
  }

  grid.reset_molecule_tile(sm.s.grid_tile_index);
  new_grid.set_molecule_tile(new_tile_index, sm.id);
  sm.s.grid_tile_index = new_tile_index;
  sm.s.wall_index = new_wall_index;

  sm.s.pos = new_loc;

  return true;
}


static void tiny_diffuse_3D(
    Partition& p,
    Molecule& vm,
    const Vec3& displacement,
    const wall_index_t previous_reflected_wall,
    Vec3& new_pos) {

  assert(vm.is_vol());
  assert(vm.v.subpart_index != SUBPART_INDEX_INVALID && "Molecule must be already placed into a subvolume");

  const BNG::Species& species = p.get_species(vm.species_id);
  Vec3 temp_displacement = displacement;
  subpart_index_t new_subpart_index;
  CollisionsVector collisions;

  // remember original location
  new_pos = vm.v.pos;

  // NOTE: can be optimized by ignoring molecule collisions
  // changes vm.pos
  ray_trace_vol(
        p, p.aux_rng,
        vm.id, species.can_vol_react(),
        previous_reflected_wall,
        temp_displacement,
        collisions
  );
  assert(vm.v.subpart_index == p.get_subpart_index(vm.v.pos));

  // sort collisions by time
  sort_collisions_by_time(collisions);

  Vec3 new_displacement = displacement;
  for (size_t collision_index = 0; collision_index < collisions.size(); collision_index++) {
    Collision& collision = collisions[collision_index];

    // stop after first collision
    if (collision.is_wall_collision() && cmp_le(collision.time, 1, EPS)) {
      // different implementation than in MCell3
      new_displacement = displacement * Vec3(collision.time*0.5);
    }
  }

  new_pos = new_pos + new_displacement;
}


static void reflect_absorb_check_wall(
    Partition& p,
    const Molecule& sm,
    const Wall& wall,
    bool& reflect_now,
    BNG::RxnClass*& absorb_now_rxn_class // set no non-null value when molecule should be absorbed
) {
  absorb_now_rxn_class = nullptr;

  BNG::RxnClassesVector matching_rxns;
  RxnUtils::trigger_intersect(
      p, sm, sm.s.orientation, wall, true,
      matching_rxns
  );

  // check if this wall has any reflective or absorptive region borders for
  // this molecule (aka special reactions)
  for (BNG::RxnClass* rxn_class: matching_rxns) {
    if (rxn_class->is_reflect_type()) {
      // check for REFLECTIVE border
      reflect_now = true;
      break;
    }
    else if (rxn_class->is_absorb_region_border_type_incl_all_molecules()) {
      // check for ABSORPTIVE border including a special case ALL_MOLECULES + surf class -> 0
      absorb_now_rxn_class = rxn_class;
      break;
    }
  }
  /* count hits if we absorb or reflect */
}


/*************************************************************************
reflect_absorb_inside_out:
  In: world: simulation state
      sm: molecule that is moving
      hd_head: region border hit data information
      rx: the type of reaction if any - absorptive/reflective
      matching_rxns: an array of possible reactions
      boundary_pos: the uv coordinates where we hit
      this_wall: the wall that we are on
      index_edge_was_hit: the index of the edge we just hit (0,1,2)
      reflect_now: should the sm reflect
      absorb_now: should the sm be absorbed
      this_wall_edge_region_border:
  Out: 1 if we are about to reflect or absorb (reflect_now, absorb_now). 0
       otherwise. hd_head and this_wall_edge_region_border are updated.
*************************************************************************/
static void reflect_absorb_inside_out(
    Partition& p,
    const Molecule& sm,
    const Wall& this_wall,
    const edge_index_t edge_index_that_was_hit,
    bool& reflect_now,
    BNG::RxnClass*& absorb_now_rxn_class
) {
  absorb_now_rxn_class = nullptr;

  // missing hit from the second side
  if (WallUtils::is_wall_edge_region_border(p, this_wall, edge_index_that_was_hit, true)) {
    reflect_absorb_check_wall(p, sm, this_wall, reflect_now, absorb_now_rxn_class);
  }
}


/*************************************************************************
reflect_absorb_outside_in:
  In: world: simulation state
      sm: molecule that is moving
      hd_head: region border hit data information
      rx: the type of reaction if any - absorptive/reflective
      matching_rxns: an array of possible reactions
      boundary_pos: the uv coordinates where we hit
      target_wall: the wall we hit
      this_wall: the wall that we are on
      reflect_now: should the sm reflect
      absorb_now: should the sm be absorbed
      this_wall_edge_region_border:
  Out: 1 if we are about to reflect or absorb (reflect_now, absorb_now). 0
       otherwise. hd_head is updated.
*************************************************************************/
static void reflect_absorb_outside_in(
    Partition& p,
    const Molecule& sm,
    const Wall& target_wall,
    const Wall& this_wall,
    bool& reflect_now,
    BNG::RxnClass*& absorb_now_rxn_class) {
  absorb_now_rxn_class = nullptr;

  /* index of the shared edge in the coordinate system of target wall */
  edge_index_t target_edge_index = WallUtils::find_shared_edge_index_of_neighbor_wall(this_wall, target_wall);

  if (WallUtils::is_wall_edge_region_border(p, target_wall, target_edge_index, true)) {
    reflect_absorb_check_wall(p, sm, target_wall, reflect_now, absorb_now_rxn_class);
  }
}

} // namespace DiffusionUtil

} // namespace mcell

