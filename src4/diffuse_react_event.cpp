/******************************************************************************
 *
 * Copyright (C) 2019 by
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
// TODO_LATER: optimization for diffusion constant 0, e.g. test 1172


#include <iostream>
#include <sstream>
#include <algorithm>
#include <boost/container/flat_set.hpp>

//extern "C" {
#include "rng.h"
#include "mcell_structs.h"
#include "logging.h"
//}

#include "diffuse_react_event.h"
#include "defines.h"
#include "world.h"
#include "partition.h"
#include "geometry.h"
#include "debug_config.h"
#include "debug.h"

// include implementations of utility functions
#include "geometry_utils.inc"
#include "reaction_utils.inc"
#include "collision_utils.inc"
#include "exact_disk_utils.inc"
#include "diffusion_utils.inc"
#include "grid_utils.inc"
#include "wall_utils.inc"

using namespace std;

namespace MCell {

void DiffuseReactEvent::step() {
  assert(event_time == world->get_current_iteration() && "In a normal case, event time should correspond to the iteration");
  assert(world->get_partitions().size() == 1 && "Must extend cache to handle multiple partitions");

  // for each partition
  for (Partition& p: world->get_partitions()) {
    // diffuse molecules from volume_molecule_indices_per_time_step that have the current diffusion_time_step
    uint time_step_index = p.get_molecule_list_index_for_time_step(diffusion_time_step);
    if (time_step_index != TIME_STEP_INDEX_INVALID) {
      diffuse_molecules(p, p.get_molecule_ids_for_time_step_index(time_step_index));
    }
  }
}


void DiffuseReactEvent::diffuse_molecules(Partition& p, const std::vector<molecule_id_t>& molecule_ids) {

  // we need to strictly follow the ordering in mcell3, therefore steps 2) and 3) do not use the time
  // for which they were scheduled but rather simply the order in which these "microevents" were created

  // 1) first diffuse already existing molecules
  uint existing_mols_count = molecule_ids.size();
  for (uint i = 0; i < existing_mols_count; i++) {
    molecule_id_t id = molecule_ids[i];
    // existing molecules - simulate whole time step
    diffuse_single_molecule(p, id, event_time, WallTileIndexPair());
  }


  // 2) simulate remaining time of molecules created with reactions or
  // scheduled unimolecular reactions
  // need to call .size() each iteration because the size can increase,
  // again, we are using it as a queue and we do not follow the time when
  // they were created
  for (uint i = 0; i < new_diffuse_or_unimol_react_actions.size(); i++) {
    const DiffuseOrUnimolRxnAction& action = new_diffuse_or_unimol_react_actions[i];

    assert(action.scheduled_time >= event_time && action.scheduled_time <= event_time + diffusion_time_step);

    if (action.type == DiffuseOrUnimolRxnAction::Type::DIFFUSE) {
      diffuse_single_molecule(
          p, action.id,
          action.scheduled_time,
          action.where_created_this_iteration // making a copy of this pair
      );
    }
    else {
      bool diffuse_right_away = react_unimol_single_molecule(p, action.id, action.scheduled_time, action.unimol_rx);

      // get a fresh action reference - unimol reaction could have added some items into the array
      const DiffuseOrUnimolRxnAction& new_ref_action = new_diffuse_or_unimol_react_actions[i];
      if (diffuse_right_away) {
        // if the molecule survived (e.g. in rxn like A -> A + B), then
        // diffuse it right away
        // (another option is to put it into the new_diffuse_or_unimol_react_actions
        //  but MCell3 diffuses them right away)
        diffuse_single_molecule(
            p, new_ref_action.id,
            new_ref_action.scheduled_time,
            new_ref_action.where_created_this_iteration // making a copy of this pair
        );
      }
    }
  }

  new_diffuse_or_unimol_react_actions.clear();
}


void DiffuseReactEvent::diffuse_single_molecule(
    Partition& p,
    const molecule_id_t m_id,
    const float_t diffusion_start_time, // time for which was this event scheduled
    WallTileIndexPair wall_tile_pair_where_created_this_iteration // set only for newly created molecules
) {
  assert(diffusion_start_time < event_time + diffusion_time_step);

  Molecule& m = p.get_m(m_id);

  if (m.is_defunct())
    return;


  float_t time_up_to_event_end =  event_time + diffusion_time_step - diffusion_start_time;

  // if the molecule is a "newbie", its unimolecular reaction was not yet scheduled
  if (m.has_flag(MOLECULE_FLAG_RESCHEDULE_UNIMOL_RX) != 0) {
    create_unimol_rx_action(p, m, time_up_to_event_end);
    m.clear_flag(MOLECULE_FLAG_RESCHEDULE_UNIMOL_RX);
  }

  // schedule unimol action if it is supposed to be executed in this timestep
  assert(m.unimol_rx_time == TIME_INVALID || m.unimol_rx_time >= event_time);
  if (m.unimol_rx_time != TIME_INVALID && m.unimol_rx != nullptr && m.unimol_rx_time < event_time + diffusion_time_step) {

    // now, there are two queues - local for this timestep
    // and global in partition for the following timesteps
    DiffuseOrUnimolRxnAction unimol_react_action(
        DiffuseOrUnimolRxnAction::Type::UNIMOL_REACT, m.id, m.unimol_rx_time, m.unimol_rx);
    // handle this iteration
    new_diffuse_or_unimol_react_actions.push_back(unimol_react_action);
  }
  /*else {
    p.add_unimolecular_action(diffusion_time_step, unimol_react_action);
  }*/


#ifdef DEBUG_DIFFUSION
  const Species& debug_species = p.all_species.get(m.species_id);
  float_t event_time_end = event_time + diffusion_time_step;
  DUMP_CONDITION4(
    // the subtraction of diffusion_time_step doesn't make much sense but is needed to make the dump the same as in mcell3
    // need to check it further
    const char* title =
        (debug_species.can_diffuse()) ?
            (m.is_vol() ? "Diffusing vm:" : "Diffusing sm:") :
            (m.is_vol() ? "Not diffusing vm:" : "Not diffusing sm:");
    m.dump(p, "", title, world->get_current_iteration(), diffusion_start_time);
  );
#endif

  // max_time is the time for which we should simulate the diffusion
  float_t max_time = time_up_to_event_end;
  if (m.unimol_rx_time != TIME_INVALID && m.unimol_rx_time < diffusion_start_time + max_time) {
    assert(m.unimol_rx_time >= diffusion_start_time);
    max_time = m.unimol_rx_time - diffusion_start_time;
  }


  // we might need to adjust remaining time step if this molecule has a unimolecular reaction
  // within this event's time step range
  /*float_t remaining_time_step;
  if (m.unimol_rx_time != TIME_INVALID && m.unimol_rx_time < event_time_end) { // unlikely
    assert(m.unimol_rx_time >= event_time && "Missed unimol rx");

    // the value of remaining_time_step passed as argument is time up to event_time_end
    float_t prev_time_from_event_start = diffusion_time_step - time_up_to_event_end;
    float_t new_time_from_event_start = m.unimol_rx_time - event_time;
    assert(new_time_from_event_start >= prev_time_from_event_start && "Unimol rx cannot be scheduled to the past");

    remaining_time_step = new_time_from_event_start - prev_time_from_event_start;
  }
  else {
    remaining_time_step = time_up_to_event_end;
  }
  assert(remaining_time_step <= diffusion_time_step && "remaining_time_step is how much time we should simulate in this step");
   */



  if (m.is_vol()) {
    diffuse_vol_molecule(
        p, m_id,
        max_time,
        diffusion_start_time,
        wall_tile_pair_where_created_this_iteration
    );
  }
  else {
    diffuse_surf_molecule(
        p, m_id, max_time, diffusion_start_time
    );
  }
}

// ---------------------------------- volume diffusion ----------------------------------

void sort_collisions_by_time(collision_vector_t& molecule_collisions) {
  sort( molecule_collisions.begin(), molecule_collisions.end(),
      [ ]( const auto& lhs, const auto& rhs )
      {
        assert((lhs.type != CollisionType::VOLMOL_SURFMOL && lhs.type != CollisionType::SURFMOL_SURFMOL) &&
            "Ray trace can return only vol-wall or vol-vol collisions");

        if (lhs.time < rhs.time) {
          return true;
        }
        else if (lhs.time > rhs.time) {
          return false;
        }
        else if (lhs.type == CollisionType::VOLMOL_VOLMOL && rhs.type == CollisionType::VOLMOL_VOLMOL) {
          // mcell3 returns collisions with molecules ordered descending by the molecule index
          // we need to maintain this behavior (needed only for identical results)
          return lhs.colliding_molecule_id > rhs.colliding_molecule_id;
        }
        else {
          return false;
        }
      }
  );
}

void DiffuseReactEvent::diffuse_vol_molecule(
    Partition& p,
    const molecule_id_t vm_id,
    const float_t max_time,
    const float_t diffusion_start_time,
    WallTileIndexPair& wall_tile_pair_where_created_this_iteration
) {
  Molecule& m = p.get_m(vm_id);
  const Species& species = p.all_species.get(m.species_id);

  // diffuse each molecule - get information on position change
  vec3_t displacement;

  float_t steps = 1.0;
  float_t t_steps = 1.0;
  float_t rate_factor = 1.0;
  float_t r_rate_factor = 1.0;
  DiffusionUtil::compute_vol_displacement(
      species, max_time, world->rng,
      displacement, rate_factor, r_rate_factor, steps, t_steps
  );

#ifdef DEBUG_TIMING
  DUMP_CONDITION4(
      dump_vol_mol_timing(
          "- Timing vm", p.stats.get_current_iteration(), vm_id,
          diffusion_start_time, max_time, m.unimol_rx_time,
          rate_factor, r_rate_factor, steps, t_steps
      );
  );
#endif

#ifdef DEBUG_DIFFUSION
  DUMP_CONDITION4(
      if (species.can_diffuse()) {
        displacement.dump("  displacement:", "");
      }
  );
#endif
  // note: we are ignoring use_expanded_list setting compared to mcell3

  // cut the displacement so that it does not cross partition?
  // or how to simplify this case when we know that we will definitely
  // be stopped by a wall (we need to fail if not anyway)

  vec3_t remaining_displacement = displacement;

  RayTraceState state;
  collision_vector_t molecule_collisions;
  bool was_defunct = false;
  vec3_t new_pos;
  subpart_index_t new_subpart_index;
  wall_index_t reflected_wall_index = WALL_INDEX_INVALID;

  //float_t updated_remaining_time_step = remaining_time_step; // == t_steps

  // ????
  float_t elapsed_molecule_time = diffusion_start_time; // == vm->t

  do {
    state =
        ray_trace_vol(
            p, world->rng,
            m /* changes position */,
            reflected_wall_index,
            remaining_displacement,
            molecule_collisions,
            new_pos,
            new_subpart_index
        );

    sort_collisions_by_time(molecule_collisions);

#ifdef DEBUG_COLLISIONS
    DUMP_CONDITION4(
        Collision::dump_array(p, molecule_collisions);
    );
#endif

    // evaluate and possible execute collisions and reactions
    for (size_t collision_index = 0; collision_index < molecule_collisions.size(); collision_index++) {
      Collision& collision = molecule_collisions[collision_index];

      assert(collision.time >= 0 && collision.time <= 1);

      if (collision.is_vol_mol_vol_mol_collision()) {
        // ignoring immediate collisions
        if (CollisionUtil::is_immediate_collision(collision.time)) {
          continue;
        }

        // evaluate reaction associated with this collision
        // for now, do the change right away, but we might need to cache these changes and
        // do them after all diffusions were finished
        // warning: might invalidate references to p.volume_molecules array! returns true in that case
        // also, if this is a reaction where this diffused product is kept, we simply continue with diffusion
        if (collide_and_react_with_vol_mol(
              p, collision, remaining_displacement,
              t_steps, r_rate_factor, elapsed_molecule_time)
        ) {
          // molecule was destroyed
          was_defunct = true;
          break;
        }
      }
      else if (collision.is_wall_collision()) {
        Molecule& vm_new_ref = p.get_m(vm_id);

        Wall& colliding_wall = p.get_wall(collision.colliding_wall_index);

        // call callback if the user registered one
        wall_hit_callback_func wall_hit_callback = world->get_wall_hit_callback();
        if (wall_hit_callback != nullptr) {
          WallHitInfo info;
          info.molecule_id = vm_new_ref.id;
          info.geometry_object_id = colliding_wall.object_id;
          info.wall_id = colliding_wall.id;
          info.time = event_time + collision.time;
          info.pos = collision.pos;

          wall_hit_callback(info, world->wall_hit_callback_clientdata);
        }
#ifdef DEBUG_WALL_COLLISIONS
        cout << "Wall collision: \n";
        const GeometryObject* geom_obj = p.get_geometry_object_if_exists(colliding_wall.object_id);
        assert(geom_obj != nullptr);
        cout << "  mol id: " << vm_new_ref.id << ", species: " << p.all_species.get(vm_new_ref.species_id).name << "\n";
        cout << "  obj: " << geom_obj->name << ", id: " << geom_obj->id << "\n";
        cout << "  wall id: " << colliding_wall.id << "\n";
        cout << "  time: " << float_t(world->get_current_iteration()) + collision.time << "\n";
        cout << "  pos: " << collision.pos << "\n";
#endif

        // check possible reaction with surface molecules
        if (species.has_flag(SPECIES_FLAG_CAN_VOLSURF) && colliding_wall.has_initialized_grid()) {
          int collide_res = collide_and_react_with_surf_mol(
              p, collision, t_steps,
              r_rate_factor, elapsed_molecule_time,
              wall_tile_pair_where_created_this_iteration
          );

          if (collide_res == 1) { // FIXME: use enum
            was_defunct = true;
            return;
          }
        }

        if (!was_defunct) {
          elapsed_molecule_time += t_steps * collision.time;
          int res = CollisionUtil::reflect_or_periodic_bc(
              p, collision,
              vm_new_ref, remaining_displacement, t_steps, reflected_wall_index
          );
          assert(res == 0 && "Periodic box BCs are not supported yet");
        }

        // molecule could have been moved
        subpart_index_t subpart_after_wall_hit = p.get_subpart_index(vm_new_ref.v.pos);
        // change subpartition if needed
        p.change_molecule_subpartition(vm_new_ref, subpart_after_wall_hit);

        break; // we reflected and did not react, do ray_trace again
      }
    }

  } while (unlikely(state != RayTraceState::FINISHED && !was_defunct));

  if (!was_defunct) { // FIXME: unify - we are changing the partition, but not the wall here - probably move to diffuse diffuse_vol_molecule
    // need to get a new reference
    Molecule& m_new_ref = p.get_m(vm_id);

    if (m_new_ref.is_vol()) {
      // finally move molecule to its destination
      m_new_ref.v.pos = new_pos;

#ifdef DEBUG_DIFFUSION
      DUMP_CONDITION4(
        // the subtraction of diffusion_time_step doesn't make much sense but is needed to make the dump the same as in mcell3
        // need to check it further
        m.dump(p, "", "- Final vm position:", world->get_current_iteration(), /*event_time_end*/ 0);
      );
#endif

      // are we still in the same partition or do we need to move?
      bool move_to_another_partition = !p.in_this_partition(m_new_ref.v.pos);
      if (move_to_another_partition) {
        mcell_log("Error: Crossing partitions is not supported yet.\n");
        exit(1);
      }

      // change subpartition
      p.change_molecule_subpartition(m_new_ref, new_subpart_index);
    }
  }
}

// move to collision utils or to partition
static vec3_t get_displacement_up_to_partition_boundary(
    const Partition& p,
    const vec3_t& pos,
    const vec3_t& displacement
    ) {

  vec3_t new_pos = pos + displacement;
  assert(p.in_this_partition(pos));
  assert(!p.in_this_partition(new_pos));

  // TODO: similar code is in collect_crossed_subparts, clean it up
  // direction of the displacement
  // TODO: get rid of the zero div check??

  vec3_t displacement_nonzero = displacement;
  debug_guard_zero_div(displacement_nonzero);
  ivec3_t dir_urb_direction = ivec3_t(glm::greaterThan(displacement_nonzero, vec3_t(0)));

  ivec3_t dir_urb_addend;
  dir_urb_addend.x = (dir_urb_direction.x == 0) ? -1 : 1;
  dir_urb_addend.y = (dir_urb_direction.y == 0) ? -1 : 1;
  dir_urb_addend.z = (dir_urb_direction.z == 0) ? -1 : 1;

  // position of edges in our direction
  vec3_t partition_edges =
      p.get_origin_corner()
      + vec3_t(dir_urb_direction) * p.config.partition_edge_length;

  vec3_t diff = partition_edges - pos;

  vec3_t displacement_rcp = 1.0/displacement_nonzero; // POSSIBLE ZERO DIV

  // time we hit a boundary
  float_t hit_time = 1;



  // first check whether we are not in fact touching one of the boundaries
  if (abs(diff.x) < EPS) {
    // only update the xyz subpartition index
    // curr_subpart_indices.x += dir_urb_addend.x;
    assert(false && "TODO");
  }
  else if (abs(diff.y) < EPS) {
    //curr_subpart_indices.y += dir_urb_addend.y;
    assert(false && "TODO");
  }
  else if (abs(diff.z) < EPS) {
    //curr_subpart_indices.z += dir_urb_addend.z;
    assert(false && "TODO");
  }
  else {
    // compute time for the next subpartition collision, let's assume that displacemnt
    // is our speed vector and the total time to travel is 1
    //
    // pos(time) = pos + displacement * time, therefore
    // time = (pos(time) - vm.v.pos) / displacement
    // =>
    // time_to_subpart_edge = (subpart_edge - vm.v.pos) / displacement_speed
    vec3_t coll_times = diff * displacement_rcp;
    assert(coll_times.x >= 0 && coll_times.y >= 0 && coll_times.z >= 0
      && "Subpartition 'edges' must be computed from direction, we cannot hit a subpart boundary that is behind us");

    // which of the times is the smallest? - i.e. which boundary we hit first
    if (coll_times.x >= 0 && coll_times.x < coll_times.y && coll_times.x <= coll_times.z) {
      // x
      hit_time = coll_times.x;
    }
    else if (coll_times.y >= 0 && coll_times.y <= coll_times.z) {
      // y
      hit_time = coll_times.y;
    }
    else if (coll_times.z >= 0) {
      // z
      hit_time = coll_times.z;
    }
    else {
      assert(false && "Collision time must not be negative");
    }
  }

  // there might be some floating point imprecisions, we want this value to be clearly in our partition,
  // so let's make the time a bit smaller
  // the displacement value is used only to find out which subpartitions we are crossing
  vec3_t new_displacement = displacement * (hit_time - EPS);
  assert(p.in_this_partition(pos + new_displacement));

  return new_displacement;
}


// collect possible collisions for molecule vm that has to displace by remaining_displacement,
// returns possible collisions in molecule_collisions, new position in new_pos and
// index of the new subparition in new_subpart_index
// later, this will check collisions until a wall is hit
// we assume that wall collisions do not occur so often
RayTraceState ray_trace_vol(
    Partition& p,
    rng_state& rng,
    Molecule& vm, // molecule that we are diffusing, we are changing its pos  and possibly also subvolume
    const wall_index_t previous_reflected_wall, // is WALL_INDEX_INVALID when our molecule did not replect from anything this iddfusion step yet
    vec3_t& remaining_displacement, // in/out - recomputed if there was a reflection
    collision_vector_t& collisions, // both mol mol and wall collisions
    vec3_t& new_pos,
    subpart_index_t& new_subpart_index
    ) {
  p.stats.inc_ray_voxel_tests();

  RayTraceState res_state = RayTraceState::FINISHED;
  collisions.clear();

  float_t radius = p.config.rx_radius_3d;

  // if we would get out of this partition, cut off the displacement
  // so we check collisions just here
  vec3_t partition_displacement;
  if (!p.in_this_partition(vm.v.pos + remaining_displacement)) {
    partition_displacement = get_displacement_up_to_partition_boundary(p, vm.v.pos, remaining_displacement);
  }
  else {
    partition_displacement = remaining_displacement;
  }

  // first get what subpartitions might be relevant
  SubpartIndicesVector crossed_subparts_for_walls;
  subpart_indices_set_t crossed_subparts_for_molecules;
  subpart_index_t last_subpartition_index;
  CollisionUtil::collect_crossed_subparts(
      p, vm, partition_displacement,
      radius, p.config.subpartition_edge_length,
      true, crossed_subparts_for_walls,
      crossed_subparts_for_molecules, last_subpartition_index
  );


  // changed when wall was hit
  vec3_t displacement_up_to_wall_collision = remaining_displacement;
  vec3_t& corrected_displacement = remaining_displacement;

  // check wall collisions in the crossed subparitions,
  if (!crossed_subparts_for_walls.empty()) {
    for (subpart_index_t subpart_index: crossed_subparts_for_walls) {

      CollisionUtil::collect_wall_collisions( // mcell3 does this only for the current subvol
          p,
          vm,
          subpart_index,
          previous_reflected_wall,
          rng,
          corrected_displacement,
          displacement_up_to_wall_collision, // may be update in case we need to 'redo' the collision detection
          collisions
      );

      // stop at first crossing because crossed_subparts_for_walls are ordered
      // and we are sure that if we hit a wall in the actual supartition, we cannot
      // possibly hit another wall in a subparition that follows
      if (!collisions.empty()) {
        res_state = RayTraceState::RAY_TRACE_HIT_WALL;
        break;
      }
    }
  }

  if (res_state == RayTraceState::RAY_TRACE_HIT_WALL) {
    // recompute collect_crossed_subparts if there was a wall collision
    // NOTE: this can be in theory done more efficiently if we knew the order of subpartitions that we hit in the previous call
    crossed_subparts_for_molecules.clear();
    CollisionUtil::collect_crossed_subparts(
        p, vm, displacement_up_to_wall_collision,
        radius,
        p.config.subpartition_edge_length,
        false, crossed_subparts_for_walls, // not filled this time
        crossed_subparts_for_molecules, last_subpartition_index
    );
  }

  // check molecule collisions for each SP
  for (subpart_index_t subpart_index: crossed_subparts_for_molecules) {
    // get cached reacting molecules for this SP
    uint_set<molecule_id_t>& sp_reactants = p.get_volume_molecule_reactants(subpart_index, vm.species_id);

    // for each molecule in this SP
    for (molecule_id_t colliding_vm_id: sp_reactants) {
      CollisionUtil::collide_mol_loop_body(
          p,
          vm,
          colliding_vm_id,
          corrected_displacement,// needs the full displacement to compute reaction time displacement_up_to_wall_collision,
          radius,
          collisions
      );
    }
  }

  // these values are valid only when RAY_TRACE_FINISHED is returned
  new_subpart_index = last_subpartition_index;
  new_pos = vm.v.pos + remaining_displacement;

  return res_state; // no wall was hit
}


// handle collision of two volume molecules: checks probability of reaction,
// executes this reaction, removes reactants and creates products
// returns true if the first reactant was destroyed
bool DiffuseReactEvent::collide_and_react_with_vol_mol(
    Partition& p,
    const Collision& collision,
    vec3_t& displacement,
    const float_t t_steps,
    const float_t r_rate_factor,
    const float_t elapsed_molecule_time
)  {

  Molecule& colliding_molecule = p.get_m(collision.colliding_molecule_id); // am
  Molecule& diffused_molecule = p.get_m(collision.diffused_molecule_id); // m

  // returns 1 when there are no walls at all
  float_t factor = ExactDiskUtil::exact_disk(
      p, collision.pos, displacement, p.config.rx_radius_3d,
      diffused_molecule, colliding_molecule,
      p.config.use_expanded_list
  );

  if (factor < 0) { /* Probably hit a wall, might have run out of memory */
    return 0; /* Reaction blocked by a wall */
  }


  // TODO: unify usage to rxn and rxn class pointers,
  //       they might not be set in Collision struct, so pointers are the only way
  const RxnClass* rxn_class = collision.rxn_class;
  assert(rxn_class != nullptr);

  //  rx->prob_t is always NULL in out case update_probs(world, rx, m->t);
  // returns which reaction pathway to take
  float_t scaling = factor * r_rate_factor;
  int i = RxUtil::test_bimolecular(
      rxn_class, world->rng, colliding_molecule, diffused_molecule, scaling, 0);

  if (i < RX_LEAST_VALID_PATHWAY) {
    return false;
  }
  else {
    // might invalidate references
    int j = outcome_bimolecular(p, collision, i, elapsed_molecule_time + t_steps * collision.time);
    assert(j == RX_DESTROY || j == RX_A_OK);
    return j == RX_DESTROY;
  }
}

/******************************************************************************
 *
 * collide_and_react_with_surf_mol is a helper function used in diffuse_3D to
 * handle collision of a diffusing 3D molecule with a surface molecule
 *
 * Return values:
 *
 * -1 : nothing happened - continue on with next smash targets
 *  0 : reaction happened and we still exist but are done with the current smash
 *  1 : reaction happened and we are destroyed
 *      target
 *
 ******************************************************************************/
int DiffuseReactEvent::collide_and_react_with_surf_mol(
    Partition& p,
    const Collision& collision,
    const float_t remaining_time_step,
    const float_t r_rate_factor,
    const float_t elapsed_molecule_time,
    WallTileIndexPair& where_created_this_iteration
) {
  Wall& wall = p.get_wall(collision.colliding_wall_index);
  Grid& grid = wall.grid;

  tile_index_t j = GridUtil::xyz2grid_tile_index(p, collision.pos, wall);
  assert(j != TILE_INDEX_INVALID);

  molecule_id_t colliging_mol_id = grid.get_molecule_on_tile(j);
  if (colliging_mol_id == MOLECULE_ID_INVALID) {
    return -1;
  }

  orientation_t collision_orientation = ORIENTATION_DOWN;
  if (collision.type == CollisionType::WALL_FRONT) {
    collision_orientation = ORIENTATION_UP;
  }

  Molecule& colliding_molecule = p.get_m(colliging_mol_id);
  assert(colliding_molecule.is_surf());

  Molecule& diffused_molecule = p.get_m(collision.diffused_molecule_id); // m
  assert(diffused_molecule.is_vol());

  // Avoid rebinding - i.e. when we created a volume molecule from a surf+vol->surf+vol
  // reaction, we do not want the molecule to react again
  if (where_created_this_iteration.wall_index == wall.index &&
      where_created_this_iteration.tile_index == j) {
    // However, let this occur next time
    where_created_this_iteration.wall_index = WALL_INDEX_INVALID;
    where_created_this_iteration.tile_index = TILE_INDEX_INVALID;
    return -1;
  }

  RxnClassesVector matching_rxn_classes;
  RxUtil::trigger_bimolecular(
    p.all_reactions.bimolecular_reactions_map,
    diffused_molecule, colliding_molecule,
    collision_orientation, colliding_molecule.s.orientation,
    matching_rxn_classes
  );

  if (matching_rxn_classes.empty()) {
    return -1;
  }

  assert(matching_rxn_classes.size() == 1 && "There should be max 1 rxn class");

  // FIXME: this code is very similar to code in react_2D_all_neighbors
  small_vector<float_t> scaling_coefs;
  for (size_t i = 0; i < matching_rxn_classes.size(); i++) {
    const RxnClass* rxn = matching_rxn_classes[i];
    assert(rxn != nullptr);

    scaling_coefs.push_back(r_rate_factor / grid.binding_factor);
  }

  int selected_rx_pathway;
  int reactant_index;
  if (matching_rxn_classes.size() == 1) {
    selected_rx_pathway = RxUtil::test_bimolecular(
        matching_rxn_classes[0], world->rng,
        diffused_molecule, colliding_molecule,
        scaling_coefs[0], 0);

    assert(selected_rx_pathway <= 0 && "Only one pathway supported for now (with index 0)");
    reactant_index = 0;
  }
  else {
    // TODO: cleanup
    // TODO: is test_many_bimolecular really neeeded?
    assert(false && "This should not happen - why would we need multiple rxn classes?");
    /*bool all_neighbors_flag = true;
    reactant_index = RxUtil::test_many_bimolecular(matching_rxn_classes, scaling_coefs, 0, world->rng, false);
    selected_rx_pathway = 0; // TODO_PATHWAYS: use value from test_many_bimolecular*/
  }

  // TODO: cleanup
  if (reactant_index == RX_NO_RX || selected_rx_pathway < RX_LEAST_VALID_PATHWAY) {
    return -1; /* No reaction */
  }

  assert(selected_rx_pathway == 0 && "TODO");

  /* run the reaction */
  float_t collision_time = elapsed_molecule_time + remaining_time_step * collision.time;

  Collision rx_collision = Collision(
      CollisionType::VOLMOL_SURFMOL,
      &p,
      collision.diffused_molecule_id,
      collision_time, // unused? FIXME: find places where the collision time is not used and remove
      collision.pos,
      colliding_molecule.id,
      matching_rxn_classes[0]
  );

  int outcome_bimol_result = outcome_bimolecular(
      p, rx_collision, selected_rx_pathway, collision_time
  );

  if (outcome_bimol_result == RX_DESTROY) {
    return 1;
  }
  else {
    return -1;
  }

}

// ---------------------------------- surface diffusion ----------------------------------

void DiffuseReactEvent::diffuse_surf_molecule(
    Partition& p,
    const molecule_id_t sm_id,
    const float_t max_time,
    const float_t diffusion_start_time
) {
  Molecule& sm = p.get_m(sm_id);
  const Species& species = p.all_species.get(sm.species_id);

  float_t steps = 0.0;
  float_t t_steps = 0.0;
  float_t space_factor = 0.0;
  float_t elapsed_molecule_time = 0.0;

  wall_index_t original_wall_index = sm.s.wall_index;

  /* Where are we going? */
  if (species.get_time_step() > max_time) {
    t_steps = max_time;
    steps = max_time / species.get_time_step() ;
  }
  else {
    t_steps = species.get_time_step();
    steps = 1.0;
  }
  if (steps < EPS) {
    t_steps = EPS * species.get_time_step();
    steps = EPS;
  }

  if (species.get_space_step() != 0) {

    if (steps == 1.0) {
      space_factor = species.get_space_step();
    }
    else {
      space_factor = species.get_space_step() * sqrt_f(steps);
    }

#ifdef DEBUG_TIMING
  DUMP_CONDITION4(
      dump_surf_mol_timing(
          "- Timing sm", p.stats.get_current_iteration(), sm_id,
          diffusion_start_time, max_time, sm.unimol_rx_time,
          space_factor, steps, t_steps
      );
  );
#endif

    for (int find_new_position = (SURFACE_DIFFUSION_RETRIES + 1);
         find_new_position > 0; find_new_position--) {

      vec2_t displacement;
      DiffusionUtil::compute_surf_displacement(species, space_factor, world->rng, displacement);


  #ifdef DEBUG_DIFFUSION
    DUMP_CONDITION4(
        if (species.can_diffuse()) {
          displacement.dump("  displacement:", "");
        }
    );
  #endif

      assert(!species.has_flag(SPECIES_FLAG_SET_MAX_STEP_LENGTH) && "not supported yet");

      // ray_trace does the movement and all other stuff
      vec2_t new_loc;
      wall_index_t new_wall_index =
          ray_trace_surf(p, species, sm, displacement, new_loc/*, elapsed_molecule_time*/);

      // Either something ambiguous happened or we hit absorptive border
      if (new_wall_index == WALL_INDEX_INVALID) {
  #if 0
        if (kill_me == 1) {
          /* molecule hit ABSORPTIVE region border */
          if (rxp == NULL) {
            mcell_internal_error("Error in 'ray_trace_2D()' after hitting "
                                 "ABSORPTIVE region border.");
          }
          if (hd_info != NULL) {
            count_region_border_update(world, sm->properties, hd_info, sm->id);
          }
          int result = outcome_unimolecular(world, rxp, 0,
                                        (struct abstract_molecule *)sm, sm->t);
          if (result != RX_DESTROY) {
            mcell_internal_error("Molecule should disappear after hitting "
                                 "ABSORPTIVE region border.");
          }
          delete_void_list((struct void_list *)hd_info);
          hd_info = NULL;
          return NULL;
        }

        if (hd_info != NULL) {
          delete_void_list((struct void_list *)hd_info);
          hd_info = NULL;
        }
  #endif
        continue; /* Something went wrong--try again */
      }

      // After diffusing, are we still on the SAME triangle?
      if (new_wall_index == sm.s.wall_index) {
        if (DiffusionUtil::move_sm_on_same_triangle(p, sm, new_loc)) {
          continue;
        }
      }
      // After diffusing, we ended up on a NEW triangle.
      else {
        if (DiffusionUtil::move_sm_to_new_triangle(p, sm, new_loc, new_wall_index)) {
          continue;
        }
      }
      find_new_position = 0;
    }
  } // if (species.space_step != 0)


  // NOTE: what about molecules that cannot diffuse?
  bool sm_still_exists = true;
  assert(!species.has_flag(SPECIES_FLAG_CAN_SURFSURFSURF) && "Not supported");
  if (species.has_flag(SPECIES_FLAG_CAN_SURFSURF)) {
    assert(!species.has_flag(SPECIES_FLAG_CANT_INITIATE) && "Not sure what to do here");

    // the time t_steps should tell when the reaction occurred and it is quite weird because
    // it has nothing to do with the time spent diffusing
    sm_still_exists = react_2D_all_neighbors(p, sm, t_steps, diffusion_start_time, elapsed_molecule_time);
  }

  if (sm_still_exists) {
    // reactions in react_2D_all_neighbors could have invalidated the molecules array
    Molecule& new_m_ref = p.get_m(sm_id);

    // for some reason, mcell3 defines a new unimol time if the molecule has moved
    bool changed_wall = new_m_ref.s.wall_index != original_wall_index;
    bool diffusible = species.can_diffuse();
    bool can_surf_surf_react = species.has_flag(SPECIES_FLAG_CAN_SURFSURF);

    if (diffusible || can_surf_surf_react) {

      // we don't have to remove the molecule from the schedule, we can just change its unimol_rx_time,
      // this time is checked and against the scheduled time
      // mcell3 compatibility: we might change the schedule only if it is not already scheduled for this time step
      if ((!diffusible || changed_wall) &&
          new_m_ref.unimol_rx_time >= event_time + diffusion_time_step) {
        new_m_ref.unimol_rx_time = TIME_INVALID;
        new_m_ref.set_flag(MOLECULE_FLAG_RESCHEDULE_UNIMOL_RX);
      }
    }
  }

}


// returns true if molecule survived
// TODO: merge with collide_and_react_with_surf_mol
bool DiffuseReactEvent::react_2D_all_neighbors(
    Partition& p,
    Molecule& sm,
    const float_t time, // same argument as t passed in mcell3 (come up with a better name)
    const float_t diffusion_start_time, // diffusion_start_time + elapsed_molecule_time should be the time when reaction occurred
    const float_t elapsed_molecule_time
) {
  assert(elapsed_molecule_time == 0 && "This is weird - mcell3 does not care about the time of diffusion on surface when creating products");

#ifdef DEBUG_TIMING
  DUMP_CONDITION4(
      dump_react_2D_all_neighbors_timing(
          time, diffusion_start_time + elapsed_molecule_time
      );
  );
#endif

  const Wall& wall = p.get_wall(sm.s.wall_index);

  TileNeighborVector neighbors;
  GridUtil::find_neighbor_tiles(p, sm, wall, sm.s.grid_tile_index, false, /* true,*/ neighbors);

  if (neighbors.empty()) {
    return true;
  }

  const Species& sm_species = p.all_species.get(sm.species_id);


  size_t l = 0;
  // array, each item corresponds to one potential reaction
  small_vector<float_t> correction_factors;
  small_vector<molecule_id_t> reactant_molecule_ids;
  RxnClassesVector matching_rxn_classes;

  /* step through the neighbors */
  for (const WallTileIndexPair& neighbor: neighbors) {
    Wall& nwall = p.get_wall(neighbor.wall_index);
    Grid& ngrid = nwall.grid;

    // is there something on the tile?
    // TODO_LATER: this filtering should be done already while looking for neighbors
    molecule_id_t nid = ngrid.get_molecule_on_tile(neighbor.tile_index);
    if (nid == MOLECULE_ID_INVALID) {
      continue;
    }

    Molecule& nsm = p.get_m(nid);
    const Species& nsm_species = p.all_species.get(nsm.species_id);

#ifdef DEBUG_REACTIONS
    DUMP_CONDITION4(
      // the subtraction of diffusion_time_step doesn't make much sense but is needed to make the dump the same as in mcell3
      // need to check it further
      nsm.dump(p, "", "  checking in react_2D_all_neighbors: ", world->get_current_iteration(), 0.0/*event_time_end - time_up_to_event_end - diffusion_time_step*//*time ???*/);
    );
#endif

    /* check whether the neighbor molecule is behind
       the restrictive region boundary   */
    if ((sm_species.has_flag(SPECIES_FLAG_CAN_REGION_BORDER) || nsm_species.has_flag(SPECIES_FLAG_CAN_REGION_BORDER)) &&
        sm.s.wall_index != nsm.s.wall_index
    ) {
      /* INSIDE-OUT check */
      if (WallUtil::walls_belong_to_at_least_one_different_restricted_region(p, wall, sm, nwall, nsm)) {
        continue;
      }

      /* OUTSIDE-IN check */
      // note: the pairing wall is same as in mcell3, TODO: explain why is it so
      if (WallUtil::walls_belong_to_at_least_one_different_restricted_region(p, wall, nsm, nwall, sm)) {
        continue;
      }
    }

    assert(!nsm_species.has_flag(SPECIES_FLAG_EXTERNAL_SPECIES) && "TODO_LATER");

    // returns value >=1 if there can be a reaction
    size_t orig_num_rxsn = matching_rxn_classes.size();
    RxUtil::trigger_bimolecular_orientation_from_mols(
        p.all_reactions.bimolecular_reactions_map,
        sm, nsm,
        matching_rxn_classes
    );

    // extend arrays holding additional information
    // FIXME: the same code is in collide_and_react_with_surf_mol
    for (size_t i = orig_num_rxsn; i < matching_rxn_classes.size(); i++) {
      const RxnClass* rxn = matching_rxn_classes[i];
      assert(rxn != nullptr);

      correction_factors.push_back(time / ngrid.binding_factor);
      reactant_molecule_ids.push_back(nsm.id);
    }
  }

  size_t num_matching_rxns = matching_rxn_classes.size();
  if (num_matching_rxns == 0) {
    return true;
  }

  reaction_index_t selected_reaction_index;
  Collision collision;

  /* Calculate local_prob_factor for the reaction probability.
     Here we convert from 3 neighbor tiles (upper probability
     limit) to the real "num_nbrs" neighbor tiles. */
  float_t local_prob_factor = 3.0 / neighbors.size();
  int rxn_class_index;
  if (num_matching_rxns == 1) {
    // figure out what should happen
    selected_reaction_index = RxUtil::test_bimolecular(
        matching_rxn_classes[0], world->rng,
        sm, p.get_m(reactant_molecule_ids[0]),
        correction_factors[0], local_prob_factor);

    // there is just one possible class == one pair of reactants
    rxn_class_index = 0;
  }
  else {
    bool all_neighbors_flag = true;
    rxn_class_index =
        RxUtil::test_many_bimolecular(matching_rxn_classes, correction_factors, local_prob_factor, world->rng, all_neighbors_flag, selected_reaction_index);
    selected_reaction_index = 0; // TODO_PATHWAYS: use value from test_many_bimolecular
  }

  if (rxn_class_index == RX_NO_RX || selected_reaction_index < RX_LEAST_VALID_PATHWAY) {
    return true; /* No reaction */
  }

  float_t collision_time = diffusion_start_time + elapsed_molecule_time;

  collision = Collision(
      CollisionType::SURFMOL_SURFMOL,
      &p, sm.id, collision_time,
      reactant_molecule_ids[rxn_class_index],
      matching_rxn_classes[rxn_class_index]
  );

  /* run the reaction */
  int outcome_bimol_result = outcome_bimolecular(
      p, collision, selected_reaction_index, collision_time
  );

  return outcome_bimol_result != RX_DESTROY;
}






/*************************************************************************
ray_trace_2D:
  In: world: simulation state
      sm: molecule that is moving
      disp: displacement vector from current to new location
      pos: place to store new coordinate (in coord system of new wall)
      kill_me: flag that tells that molecule hits ABSORPTIVE region border
           (value = 1)
      rxp: reaction object (valid only in case of hitting ABSORPTIVE region
         border)
      hit_data_info: region border hit data information
  Out: Return wall at endpoint of movement vector. Otherwise NULL if we hit
       ambiguous location or if SM was absorbed.
       pos: location of that endpoint in the coordinate system of the new wall.
       kill_me, rxp, and hit_data_info will all be updated if we hit absorptive
       boundary.
*************************************************************************/
wall_index_t DiffuseReactEvent::ray_trace_surf(
    Partition& p,
    const Species& species,
    Molecule& sm,
    vec2_t& remaining_displacement,
    vec2_t& new_pos/*,
    float_t& elapsed_molecule_time*/
) {
  const Wall* this_wall = &p.get_wall(sm.s.wall_index);

  vec2_t orig_pos = sm.s.pos;
  vec2_t this_pos = sm.s.pos;
  vec2_t this_disp = remaining_displacement;

  /* Will break out with return or break when we're done traversing walls */
  while (1) {

    int this_wall_edge_region_border = 0;
    //bool absorb_now = 0;

    /* Index of the wall edge that the SM hits */
    vec2_t boundary_pos;
    edge_index_t edge_index_that_was_hit =
        GeometryUtil::find_edge_point(*this_wall, this_pos, this_disp, boundary_pos);

    // Ambiguous edge collision. Give up and try again from diffuse_2D.
    if (edge_index_that_was_hit == EDGE_INDEX_CANNOT_TELL) {
      sm.s.pos = orig_pos;
      // hit_data_info = hit_data_head;
      return WALL_INDEX_INVALID;
    }

    // We didn't hit the edge. Stay inside this wall. We're done!
    else if (edge_index_that_was_hit == EDGE_INDEX_WITHIN_WALL) {
      new_pos = this_pos + this_disp;

      // ???
      sm.s.pos = orig_pos;
      // *hit_data_info = hit_data_head;
      return this_wall->index;
    }


    // Neither ambiguous (EDGE_INDEX_CANNOT_TELL) nor inside wall (EDGE_INDEX_WITHIN_WALL),
    // must have hit edge (0, 1, 2)
    vec2_t old_pos = this_pos;

    /* We hit the edge - check for the reflection/absorption from the
       edges of the wall if they are region borders
       Note - here we test for potential collisions with the region
       border while moving INSIDE OUT */
    bool absorb_now = false;
    bool reflect_now = false;
    if (species.can_interact_with_border()) {
      DiffusionUtil::reflect_absorb_inside_out(
          p, sm, *this_wall, edge_index_that_was_hit,
          reflect_now, absorb_now
      );

      assert(!absorb_now && "TODO");

#if 0
      if (absorb_now) {
        assert(false && "TODO");
        /**kill_me = 1;
        *rxp = rx;
        *hit_data_info = hit_data_head;
        return NULL;*/
      }
#endif
    }

    /* no reflection - keep going */

    if (!reflect_now) {
      wall_index_t target_wall_index =
          GeometryUtil::traverse_surface(*this_wall, old_pos, edge_index_that_was_hit, this_pos);

      if (target_wall_index != WALL_INDEX_INVALID) {
        /* We hit the edge - check for the reflection/absorption from the
           edges of the wall if they are region borders
           Note - here we test for potential collisions with the region
           border while moving OUTSIDE IN */
        if (species.can_interact_with_border()) {

          const Wall& target_wall = p.get_wall(target_wall_index);
          DiffusionUtil::reflect_absorb_outside_in(
              p, sm, target_wall, *this_wall,
              reflect_now, absorb_now
          );

          assert(!absorb_now && "TODO");
        }


        if (!reflect_now) {
          this_disp = old_pos + this_disp;

          #ifndef NDEBUG
            Edge& e = const_cast<Edge&>(this_wall->edges[edge_index_that_was_hit]);
            assert(e.is_initialized());
            e.debug_check_values_are_uptodate(p);
          #endif

          vec2_t tmp_disp;
          GeometryUtil::traverse_surface(*this_wall, this_disp, edge_index_that_was_hit, tmp_disp);
          this_disp = tmp_disp - this_pos;
          this_wall = &p.get_wall(target_wall_index);
          continue;
        }
      }
    }

    /* If we reach this point, assume we reflect off the edge since there is no
     * neighboring wall
     *
     * NOTE: this_pos has been corrupted by traverse_surface; use old_pos to find
     * out whether the present wall edge is a region border
     */
    vec2_t new_disp = this_disp - (boundary_pos - old_pos);

    switch (edge_index_that_was_hit) {
      case EDGE_INDEX_0:
        new_disp.v *= -1.0;
        break;
      case EDGE_INDEX_1: {
        float_t f;
        vec2_t reflector;
        reflector.u = -this_wall->uv_vert2.v;
        reflector.v = this_wall->uv_vert2.u - this_wall->uv_vert1_u;
        f = 1.0 / sqrt_f( len2_squared(reflector) );
        reflector *= f;
        f = 2.0 * dot2(new_disp, reflector);
        new_disp -= vec2_t(f) * reflector;
        break;
      }
      case EDGE_INDEX_2: {
        float_t f;
        vec2_t reflector;
        reflector.u = this_wall->uv_vert2.v;
        reflector.v = -this_wall->uv_vert2.u;
        f = 1.0 / sqrt_f( len2_squared(reflector) );
        reflector *= f;
        f = 2.0 * dot2(new_disp, reflector);
        new_disp -= vec2_t(f) * reflector;
        break;
      }
      default:
        UNHANDLED_CASE(edge_index_that_was_hit);
    }

    this_pos.u = boundary_pos.u;
    this_pos.v = boundary_pos.v;
    this_disp.u = new_disp.u;
    this_disp.v = new_disp.v;
  } // while

  sm.s.pos = orig_pos;
}


// ---------------------------------- other ----------------------------------


void DiffuseReactEvent::create_unimol_rx_action(
    Partition& p,
    Molecule& m,
    float_t remaining_time_step
) {
  float_t curr_time = event_time + diffusion_time_step - remaining_time_step;
  assert(curr_time >= 0);

  const RxnClass* rx = RxUtil::pick_unimol_rx(world, m.species_id);
  if (rx == nullptr) {
    return;
  }

  float_t time_from_now = RxUtil::compute_unimol_lifetime(p, world->rng, m, rx);

  float_t scheduled_time = curr_time + time_from_now; // TODO: this is weird... seems mcell uses t+t2?

  // we need to store the end time to the molecule because oit is needed in diffusion to
  // figure out whether we should do the whole time step
  m.unimol_rx_time = scheduled_time;
  m.unimol_rx = rx;
}


// based on mcell3's check_for_unimolecular_reaction
// might invalidate vm references
// returns true if molecule should be diffused right away (needed for mcell3 compatibility)
bool DiffuseReactEvent::react_unimol_single_molecule(
    Partition& p,
    const molecule_id_t m_id,
    const float_t scheduled_time,
    const RxnClass* unimol_rxn_class
) {
  // the unimolecular reaction class was already selected
  assert(unimol_rxn_class != nullptr);

  Molecule& m = p.get_m(m_id);

  if (m.is_defunct()) {
    return false;
  }

  // unimolecular reactions for surface molecules can be rescheduled,
  // ignore this action in this case
  if (scheduled_time != m.unimol_rx_time) {
    return true; // TODO: check, we will probably need to return false here...
  }

  assert(scheduled_time >= event_time && scheduled_time <= event_time + diffusion_time_step);

  reaction_index_t ri = RxUtil::which_unimolecular(unimol_rxn_class, world->rng);
  return outcome_unimolecular(p, m, scheduled_time, unimol_rxn_class->get_reaction(ri));
}


// checks if reaction should probabilistically occur and if so,
// destroys reactants
// returns RX_DESTROY when the primary reactant was destroyed, RX_A_OK if the reactant A was kept
int DiffuseReactEvent::outcome_bimolecular(
    Partition& p,
    const Collision& collision,
    int path,
    float_t time // FIXME: compute time here?
) {
#ifdef DEBUG_TIMING
  DUMP_CONDITION4(
      dump_outcome_bimolecular_timing(time);
  );
#endif

  // might invalidate references

  bool keep_reacA, keep_reacB;
  int result =
      outcome_products_random(
        p,
        collision,
        time,
        path,
        keep_reacA, keep_reacB
      );

  if (result == RX_A_OK) {
    Molecule& reacA = p.get_m(collision.diffused_molecule_id);
    Molecule& reacB = p.get_m(collision.colliding_molecule_id);

#ifdef DEBUG_REACTIONS
    // reference printout first destroys B then A
    DUMP_CONDITION4(
      if (!keep_reacB) {
        reacB.dump(p, "", "  defunct m:", world->get_current_iteration(), 0, false);
      }
      if (!keep_reacA) {
        reacA.dump(p, "", "  defunct m:", world->get_current_iteration(), 0, false);
      }
    );
#endif

    // always for now
    // we used the reactants - remove them
    if (!keep_reacA) {
      p.set_molecule_as_defunct(reacA);
    }
    if (!keep_reacB) {
      p.set_molecule_as_defunct(reacB);
    }

    if (keep_reacA) {
      return RX_A_OK;
    }
    else {
      return RX_DESTROY;
    }
  }


  return result;
}

// might return RX_BLOCKED if reaction cannot occur,
// returns 0 if positions were found
int DiffuseReactEvent::find_surf_product_positions(
    Partition& p,
    const Molecule* reacA, const bool keep_reacA,
    const Molecule* reacB, const bool keep_reacB,
    const Molecule* surf_reac,
    const Rxn* rxn,
    small_vector<GridPos>& assigned_surf_product_positions) {

  uint needed_surface_positions = rxn->get_num_surf_products(p.all_species);

  small_vector<GridPos> recycled_surf_prod_positions; // this array contains information on where to place the surface products
  uint initiator_recycled_index = INDEX_INVALID;

  // find which tiles can be recycled
  if (reacA->is_surf()) {
    if (!keep_reacA) {
      recycled_surf_prod_positions.push_back( GridPos::make_with_pos(p, *reacA) );
      if (reacA->id == surf_reac->id) {
        initiator_recycled_index = 0;
      }
    }
    else if (keep_reacA) {
      // reacA is kept
      needed_surface_positions--;
    }
  }

  if (reacB != nullptr && reacB->is_surf()) {
    if (!keep_reacB) {
      recycled_surf_prod_positions.push_back( GridPos::make_with_pos(p, *reacB) );
      if (reacB->id == surf_reac->id) {
        initiator_recycled_index = recycled_surf_prod_positions.size() - 1;
      }
    }
    else if (reacB->is_surf() && keep_reacB) {
      // reacB is kept
      needed_surface_positions--;
    }
  }

  // do we need more tiles?
  TileNeighborVector vacant_neighbor_tiles;
  if (needed_surface_positions > recycled_surf_prod_positions.size()) {
    assert(surf_reac != nullptr);

    // find neighbors for the first surface reactant
    Wall& wall = p.get_wall(surf_reac->s.wall_index);
    Grid& grid = wall.grid;
    TileNeighborVector neighbor_tiles;
    GridUtil::find_neighbor_tiles(p, *surf_reac, wall, surf_reac->s.grid_tile_index, true, /* false,*/ neighbor_tiles);

    // we care only about the vacant ones (NOTE: this filtering out might be done in find_neighbor_tiles)
    // mcell3 reverses the ordering here
    for (int i = neighbor_tiles.size() - 1; i >=0; i--) {
      WallTileIndexPair& tile_info = neighbor_tiles[i];

      Grid& neighbor_grid = p.get_wall(tile_info.wall_index).grid;
      if (neighbor_grid.get_molecule_on_tile(tile_info.tile_index) == MOLECULE_ID_INVALID) {
        vacant_neighbor_tiles.push_back(tile_info);
      }
    }

    /* Can this reaction happen at all? */
    if (vacant_neighbor_tiles.size() + recycled_surf_prod_positions.size() < needed_surface_positions) {
      return RX_BLOCKED;
    }
  }

  // random assignment of positions
  uint num_tiles_to_recycle = min(rxn->products.size(), recycled_surf_prod_positions.size());
  if (num_tiles_to_recycle == 1 && recycled_surf_prod_positions.size() >= 1) {
    // NOTE: this code is overly complex and can be simplified
    if (initiator_recycled_index == INDEX_INVALID) {
      assert(recycled_surf_prod_positions.size() == 1);
      assigned_surf_product_positions[0] = recycled_surf_prod_positions[0];
    }
    else {
      assigned_surf_product_positions[0] = recycled_surf_prod_positions[initiator_recycled_index];
    }
  }
  else if (num_tiles_to_recycle > 1) {
    uint next_available_index = 0;
    uint n_players = rxn->get_num_players();
    uint n_reactants = rxn->reactants.size();

    // assign recycled positions to products
    while (next_available_index < num_tiles_to_recycle) {
      // we must have the same number of random calls as in mcell3...
      uint rnd_num = rng_uint(&world->rng) % n_players;

      // continue until we got an index of a product
      if (rnd_num < n_reactants) {
        continue;
      }

      uint product_index = rnd_num - n_reactants;
      assert(product_index < rxn->products.size());

      // we care only about surface molecules
      if (p.all_species.get(rxn->products[product_index].species_id).is_vol()) {
        continue;
      }

      // skip products that we already set
      if (assigned_surf_product_positions[product_index].initialized) {
        continue;
      }

      // set position for product with product_index
      assigned_surf_product_positions[product_index] = recycled_surf_prod_positions[next_available_index];
      next_available_index++;
    }

    // all other products are placed on one of the randomly chosen vacant tiles
    small_vector<bool> used_vacant_tiles;
    used_vacant_tiles.resize(vacant_neighbor_tiles.size(), false);

    uint n_products = rxn->products.size();
    uint num_vacant_tiles = vacant_neighbor_tiles.size();
    for (uint product_index = 0; product_index < n_products; product_index++) {

      if (assigned_surf_product_positions[product_index].initialized) {
        continue;
      }

      uint num_attempts = 0;
      bool found = false;
      while (!found && num_attempts < SURFACE_DIFFUSION_RETRIES) {

        uint rnd_num = rng_uint(&world->rng) % num_vacant_tiles;

        if (!used_vacant_tiles[rnd_num]) {
          WallTileIndexPair grid_tile_index_pair = vacant_neighbor_tiles[rnd_num];
          assigned_surf_product_positions[product_index] = GridPos::make_without_pos(p, grid_tile_index_pair);
          used_vacant_tiles[rnd_num] = true;
          found = true;
        }

        num_attempts++;
      }
      if (num_attempts >= SURFACE_DIFFUSION_RETRIES) {
        return RX_BLOCKED;
      }
    }
  }

  return 0;
}

// why is this called "random"? - check if reaction occurs is in test_bimolecular
// mcell3 version returns  cross_wall ? RX_FLIP : RX_A_OK;
// ! might invalidate references
// might return RX_BLOCKED
int DiffuseReactEvent::outcome_products_random(
    Partition& p,
    const Collision& collision,
    const float_t time,
    const reaction_index_t reaction_index,
    bool& keep_reacA,
    bool& keep_reacB
) {
#ifdef DEBUG_REACTIONS
  DUMP_CONDITION4(
      collision.dump(p, "Processing reaction:", p.stats.get_current_iteration(), event_time + time);
  );
#endif

  // TODO: unify rx vs rxn
  const Rxn* rx = nullptr;
  if (collision.is_mol_mol_reaction()) {
    const RxnClass* rxn_class = collision.rxn_class;
    assert(rxn_class != nullptr);
    rx = rxn_class->get_reaction(reaction_index);
  }
  else {
    assert(reaction_index == 0 && "For other than mol mol collision the selected reaction index must be 0");
    rx = collision.rx;
  }
  assert(rx != nullptr);
  assert(rx->reactants.size() == 1 || rx->reactants.size() == 2);

  Molecule* reacA = &p.get_m(collision.diffused_molecule_id);
  keep_reacA = false; // one product is the same as reacA
  assert(reacA != nullptr);

  Molecule* reacB = nullptr;
  keep_reacB = false; // one product is the same as reacB

  Molecule* surf_reac = nullptr;

  bool reactants_swapped = false;
  if (rx->reactants.size() == 2) {
    reacB = &p.get_m(collision.colliding_molecule_id);

    if (reacA->is_surf()) {
      surf_reac = reacA;
    }
    else if (reacB->is_surf()) {
      surf_reac = reacB;
    }

    /* Ensure that reacA and reacB are sorted in the same order as the rxn players. */
    /* Needed to maintain the same behavior as in mcell3 */
    if (SpeciesWithOrientation(reacA->species_id, reacA->get_orientation()) != rx->reactants[0]) {
      Molecule* tmp_mol = reacA;
      reacA = reacB;
      reacB = tmp_mol;
      reactants_swapped = true;
    }
    assert(rx->reactants[1].is_same_tolerate_orientation_none(reacB->species_id, reacB->get_orientation()));

    keep_reacB = rx->reactants[1].is_on_both_sides_of_rxn();
  }
  else {
    surf_reac = reacA->is_surf() ? reacA : nullptr;
  }
  assert(rx->reactants[0].is_same_tolerate_orientation_none(reacA->species_id, reacA->get_orientation()));
  keep_reacA = rx->reactants[0].is_on_both_sides_of_rxn();

  bool is_orientable = reacA->is_surf() || (reacB != nullptr && reacB->is_surf());

  WallTileIndexPair surf_reac_wall_tile;
  assert(
      ( surf_reac == nullptr ||
       (surf_reac->s.wall_index == WALL_INDEX_INVALID && surf_reac->s.grid_tile_index == TILE_INDEX_INVALID) ||
       (surf_reac->s.wall_index != WALL_INDEX_INVALID && surf_reac->s.grid_tile_index != TILE_INDEX_INVALID)) &&
      "Either both wall and tile index must be valid or both must be invalid"
  );
  if (surf_reac != nullptr) {
    // we need to remember this value now because surf_reac's grid tile is freed a bit later to make a place for
    // new product
    surf_reac_wall_tile.wall_index = surf_reac->s.wall_index;
    surf_reac_wall_tile.tile_index = surf_reac->s.grid_tile_index;
  }

  /* If the reaction involves a surface, make sure there is room for each product. */
  small_vector<GridPos> assigned_surf_product_positions; // this array contains information on where to place the surface products
  assigned_surf_product_positions.resize(rx->products.size());
  if (is_orientable) {
    int res = find_surf_product_positions(
        p, reacA, keep_reacA, reacB, keep_reacB, surf_reac, rx,
        assigned_surf_product_positions);
    if (res == RX_BLOCKED) {
      return RX_BLOCKED;
    }
  }

  // free up tiles that we are probably going to reuse
  //bool one_of_reactants_is_surf = false;
  if (reacA->is_surf() && !keep_reacA) {
    p.get_wall(reacA->s.wall_index).grid.reset_molecule_tile(reacA->s.grid_tile_index);
    reacA->s.grid_tile_index = TILE_INDEX_INVALID;
    //one_of_reactants_is_surf = true;
  }

  if (reacB != nullptr && reacB->is_surf() && !keep_reacB) {
    p.get_wall(reacB->s.wall_index).grid.reset_molecule_tile(reacB->s.grid_tile_index);
    reacB->s.grid_tile_index = TILE_INDEX_INVALID;
    //one_of_reactants_is_surf = true;
  }

  // create and place each product

  uint current_surf_product_position_index = 0;

  for (uint product_index = 0; product_index < rx->products.size(); product_index++) {
    const SpeciesWithOrientation& product = rx->products[product_index];

    // do not create anything new when the reactant is kept -
    // for bimol reactions - the diffusion simply continues
    // for unimol reactions - the unimol action action starts diffusion for the remaining timestep
    if (product.is_on_both_sides_of_rxn()) {
      // remember which reactant(s) to keep?
      continue;
    }

    const Species& species = p.all_species.get(product.species_id);

    molecule_id_t new_m_id;

    // set only for new vol mols when one of the reactants is surf, invalid by default
    WallTileIndexPair where_is_vm_created;

    float_t scheduled_time = time;
    /*if (rx->reactants.size() == 2 && species.is_vol() && !one_of_reactants_is_surf) {
      // bimolecular reaction
      // schedule new product for diffusion
      // collision.time is relative to the part that this molecule travels this diffusion step
      // so it needs to be scaled
      scheduled_time = event_time + diffusion_time_step - (time - collision.time * time);
    }
    else if (rx->reactants.size() == 2 && (species.is_surf() || one_of_reactants_is_surf)) {
      scheduled_time = event_time + collision.time;
    }
    else {
      // unimolecular reaction
      // reaction_time is the time when this new molecule was created
      scheduled_time = event_time + collision.time;
    }*/

    if (species.is_vol()) {
      // create and place a volume molecule

      Molecule vm_initialization(MOLECULE_ID_INVALID, product.species_id, collision.pos);

      // adding molecule might invalidate references of already existing molecules
      Molecule& new_vm = p.add_volume_molecule(vm_initialization);

      // id and position is used to schedule a diffusion action
      new_m_id = new_vm.id;
      if (surf_reac != nullptr) {
        where_is_vm_created = surf_reac_wall_tile;
      }

      new_vm.flags = IN_VOLUME | (species.can_diffuse() ? ACT_DIFFUSE : 0);
      new_vm.set_flag(MOLECULE_FLAG_VOL);
      new_vm.set_flag(MOLECULE_FLAG_RESCHEDULE_UNIMOL_RX);

      /* For an orientable reaction, we need to move products away from the surface
       * to ensure they end up on the correct side of the plane. */
      assert(is_orientable
          || (collision.type != CollisionType::VOLMOL_SURFMOL && collision.type != CollisionType::SURFMOL_SURFMOL)
      );
      if (is_orientable) {
        assert(surf_reac != nullptr);
        Wall& w = p.get_wall(surf_reac->s.wall_index);

        float_t bump = (product.orientation > 0) ? EPS : -EPS;
        vec3_t displacement = vec3_t(2 * bump) * w.normal;
        vec3_t new_pos_after_diffuse;

        DiffusionUtil::tiny_diffuse_3D(p, new_vm, displacement, w.index, new_pos_after_diffuse);

        // update position and subpart if needed
        new_vm.v.pos = new_pos_after_diffuse;
        subpart_index_t new_subpart = p.get_subpart_index(new_vm.v.pos);
        p.change_molecule_subpartition(new_vm, new_subpart);
      }


    #ifdef DEBUG_REACTIONS
      DUMP_CONDITION4(
        new_vm.dump(p, "", "  created vm:", world->get_current_iteration(), scheduled_time);
      );
    #endif
    }
    else {
      // see release_event_t::place_single_molecule_onto_grid, merge somehow

      // get info on where to place the product and increment the counter
      const GridPos& new_grid_pos = assigned_surf_product_positions[current_surf_product_position_index];
      current_surf_product_position_index++;
      assert(new_grid_pos.initialized);

      vec2_t pos;
      if (new_grid_pos.pos_is_set) {
        pos = new_grid_pos.pos;
      }
      else {
        // NOTE: there might be some optimization for this in mcell3
        const Wall& wall = p.get_wall(new_grid_pos.wall_index);
        pos = GridUtil::grid2uv_random(wall, new_grid_pos.tile_index, world->rng);
      }

      // create our new molecule
      Molecule sm(MOLECULE_ID_INVALID, product.species_id, pos);

      // might invalidate references of already existing molecules
      Molecule& new_sm = p.add_surface_molecule(sm);
      new_m_id = new_sm.id;
      new_sm.flags = (species.can_diffuse() ? ACT_DIFFUSE : 0) | IN_SURFACE;
      new_sm.set_flag(MOLECULE_FLAG_SURF);
      new_sm.set_flag(MOLECULE_FLAG_RESCHEDULE_UNIMOL_RX);

      // set wall and grid information
      new_sm.s.wall_index = new_grid_pos.wall_index;
      new_sm.s.grid_tile_index = new_grid_pos.tile_index;
      Grid& grid = p.get_wall(new_grid_pos.wall_index).grid;
      grid.set_molecule_tile(new_sm.s.grid_tile_index, new_sm.id);

      // and finally orientation
      new_sm.s.orientation = product.orientation;

      #ifdef DEBUG_REACTIONS
        DUMP_CONDITION4(
          new_sm.dump(p, "", "  created sm:", world->get_current_iteration(), scheduled_time);
        );
      #endif
    }

    // In this time step, we will simply simulate all results of reactions regardless on the diffusion time step of the
    // particular product
    // We always create diffuse events, unimol react events are created right before diffusion of that molecule
    // to keep MCell3 compatibility, diffusion of molecule that survived its unimol reaction must
    // be executed right away and this is handled i nDiffuseReactEvent::diffuse_molecules
    new_diffuse_or_unimol_react_actions.push_back(
        DiffuseOrUnimolRxnAction(
            DiffuseOrUnimolRxnAction::Type::DIFFUSE,
            new_m_id, scheduled_time,
            where_is_vm_created
    ));
  } // end for - product creation

  // we might need to swap info on which reactant was kept
  if (reactants_swapped) {
    bool tmp = keep_reacA;
    keep_reacA = keep_reacB;
    keep_reacB = tmp;
  }

  return RX_A_OK;
}

// ---------------------------------- unimolecular reactions ----------------------------------

// !! might invalidate references (we might reorder defuncting and outcome call later)
// returns true if molecule survived
bool DiffuseReactEvent::outcome_unimolecular(
    Partition& p,
    Molecule& m,
    const float_t scheduled_time,
    const Rxn* unimol_rx
) {
  molecule_id_t id = m.id;

  Collision collision(CollisionType::UNIMOLECULAR_VOLMOL, &p, m.id, scheduled_time, m.v.pos, unimol_rx);

  bool ignoredA, ignoredB;
  // creates new molecule(s) as output of the unimolecular reaction
  // !! might invalidate references (we might reorder defuncting and outcome call later)
  int outcome_res = outcome_products_random(p, collision, scheduled_time, 0, ignoredA, ignoredB);
  assert(outcome_res == RX_A_OK);

  Molecule& m_new_ref = p.get_m(id);

  // and defunct this molecule if it was not kept
  assert(unimol_rx->reactants.size() == 1);
  if (!unimol_rx->reactants[0].is_on_both_sides_of_rxn()) {
  #ifdef DEBUG_REACTIONS
    DUMP_CONDITION4(
      m_new_ref.dump(p, "", m_new_ref.is_vol() ? "Unimolecular vm defunct:" : "Unimolecular sm defunct:", world->get_current_iteration(), scheduled_time, false);
    );
  #endif
    p.set_molecule_as_defunct(m_new_ref);
    return false;
  }
  else {
    // we must reschedule the molecule's unimol rxn, this will happen right away 
    // during the molecule's diffusion
    m_new_ref.set_flag(MOLECULE_FLAG_RESCHEDULE_UNIMOL_RX);
    
    // a cleaner way would be to schedule a new diffusion event but we must do the diffusion right away 
    // to stay compatible with mcell3
    // float_t time_up_to_event_end = diffusion_time_step - time_from_event_start;
		// create_unimol_rx_action(p, m_new_ref, time_up_to_event_end);
    
    // molecule survived
    return true;
  }
}



// ---------------------------------- dumping methods ----------------------------------

void DiffuseReactEvent::dump(const std::string indent) {
  cout << indent << "Diffuse-react event:\n";
  std::string ind2 = indent + "  ";
  BaseEvent::dump(ind2);
  cout << ind2 << "diffusion_time_step: \t\t" << diffusion_time_step << " [float_t] \t\t\n";
}


} /* namespace mcell */
