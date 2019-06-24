/******************************************************************************
 *
 * Copyright (C) 2006-2017 by
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

#include "diffuse_react_event.h"
#include "defines.h"
#include "world.h"
#include "partition.h"
#include "geometry.h"
#include "debug_config.h"

using namespace std;

namespace mcell {

namespace diffusion_util {


// ---------------------------------- volume mol diffusion ----------------------------------

// get displacement based on scale (related to diffusion constant) and gauss random number
static void pick_vol_displacement(const float_t scale, rng_state& rng, vec3_t& displacement) {
  displacement.x = scale * rng_gauss(&rng) * .70710678118654752440;
  displacement.y = scale * rng_gauss(&rng) * .70710678118654752440;
  displacement.z = scale * rng_gauss(&rng) * .70710678118654752440;
}


// determine how far will our diffused molecule move
static void compute_vol_displacement(
    const species_t& sp,
    const float_t remaining_time_step,
    rng_state& rng,
    vec3_t& displacement,
    float_t& r_rate_factor
) {

  float_t rate_factor = (remaining_time_step == 1.0) ? 1.0 : sqrt(remaining_time_step);
  r_rate_factor = 1.0 / rate_factor;
  pick_vol_displacement(sp.space_step * rate_factor, rng, displacement);
}

// ---------------------------------- surface mol diffusion ----------------------------------

/*************************************************************************
pick_2D_displacement:
  In: v: vector2 to store the new displacement
      scale: scale factor to apply to the displacement
      rng:
  Out: No return value.  vector is set to a random orientation and a
         distance chosen from the probability distribution of a diffusing
         2D molecule, scaled by the scaling factor.
*************************************************************************/
static void pick_surf_displacement(vec2_t& v, const float_t scale, rng_state& rng) {
  const float_t one_over_2_to_16th = 1.52587890625e-5;
  vec2_t a;

  /*
   * NOTE: The below algorithm is the polar method due to Marsaglia
   * combined with a rejection method for picking uniform random
   * variates in C2.
   * Both methods are nicely described in Chapters V.4.3 and V.4.4
   * of "Non-Uniform Random Variate Generation" by Luc Devroye
   * (http://luc.devroye.org/rnbookindex.html).
   */
  float_t f;
  do {
    unsigned int n = rng_uint(&rng);

    a.u = 2.0 * one_over_2_to_16th * (n & 0xFFFF) - 1.0;
    a.v = 2.0 * one_over_2_to_16th * (n >> 16) - 1.0;
    f = len2_squared(a);
  } while ((f < EPS) || (f > 1.0));

  /*
   * NOTE: The scaling factor to go from a uniform to
   * a normal distribution is sqrt(-2log(f)/f).
   * However, since we use two normally distributed
   * variates to generate a normally distributed
   * 2d vector (with variance 1) we have to normalize
   * and divide by an additional factor of sqrt(2)
   * resulting in normalFactor.
   */
  float_t normal_factor = sqrt_f(-log_f(f) / f);
  v = a * vec2_t(normal_factor * scale);
}



/*************************************************************************
move_sm_on_same_triangle:

  This is a helper function for diffuse_2D.

  In: world: simulation state
      sm: molecule that is moving
      new_loc: this is the location we are moving to.
      previous_box: this is the periodic box we were in previously.
      new_wall: this is the new wall we ended up on
      hd_info:
  Out: The grid is created on a new triangle and we place the molecule if
       possible. Counts are updated.
*************************************************************************/
static int move_sm_on_same_triangle(
    partition_t& p,
    molecule_t& sm,
    vec2_t& new_loc
) {
  wall_t& wall = p.get_wall(sm.s.wall_index);
  grid_t& grid = p.get_grid(wall.grid_index);
  assert(grid.get_molecule_on_tile(sm.s.grid_tile_index) == sm.id);

  unsigned int new_tile_index = geom_util::uv2grid(new_loc, wall, grid);

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
      return 1; /* Pick again--full here */
    }

    grid.reset_molecule_tile(sm.s.grid_tile_index);
    grid.set_molecule_tile(new_tile_index, sm.id);
    sm.s.grid_tile_index = new_tile_index;
  }

  sm.s.pos = new_loc;
  return 0;
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
  Out: The grid is created on a new triangle and we place the molecule if
       possible. Counts are updated.
*************************************************************************/
static int move_sm_to_new_triangle(
    partition_t& p,
    molecule_t& sm,
    vec2_t& new_loc,
    const wall_index_t new_wall_index
) {
  wall_t& wall = p.get_wall(sm.s.wall_index);
  grid_t& grid = p.get_grid(wall.grid_index);
  assert(grid.get_molecule_on_tile(sm.s.grid_tile_index) == sm.id);

  wall_t& new_wall = p.get_wall(new_wall_index);

  // No SM has been here before, so we need to make a grid on this wall.
  if (new_wall.grid_index == GRID_INDEX_INVALID) {
    p.add_grid_for_wall(new_wall);
  }

  grid_t& new_grid = p.get_grid(wall.grid_index);

  /* Move to new tile */
  unsigned int new_tile_index = geom_util::uv2grid(new_loc, new_wall, new_grid);

  if (new_tile_index >= new_grid.num_tiles) {
    mcell_internal_error(
        "After ray_trace_2D to a new wall, selected u, v coordinates map "
        "to an out-of-bounds grid cell.  uv=(%.2f, %.2f) sm=%d/%d",
        new_loc.u, new_loc.v, new_tile_index, new_grid.num_tiles);
  }

  molecule_id_t molecule_id_on_tile = new_grid.get_molecule_on_tile(new_tile_index);
  if (molecule_id_on_tile != MOLECULE_ID_INVALID) {
    return 1; /* Pick again--full here */
  }

  grid.reset_molecule_tile(sm.s.grid_tile_index);
  new_grid.set_molecule_tile(new_tile_index, sm.id);
  sm.s.grid_tile_index = new_tile_index;

  sm.s.pos = new_loc;

  return 0;
}


} // namespace diffusion_util

} // namespace mcell

