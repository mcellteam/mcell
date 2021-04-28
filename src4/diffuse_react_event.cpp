/******************************************************************************
 *
 * Copyright (C) 2019,2020 by
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
// TODO: make this file shorter

#include <iostream>
#include <sstream>
#include <algorithm>
#include <boost/container/flat_set.hpp>

#include "api/mol_wall_hit_info.h"
#include "api/reaction_info.h"
#include "api/callbacks.h"

#include "rng.h"
#include "mcell_structs_shared.h"
#include "logging.h"

#include "diffuse_react_event.h"
#include "defines.h"
#include "world.h"
#include "partition.h"
#include "geometry.h"
#include "grid_position.h"
#include "region_utils.h"

#include "debug_config.h"
#include "debug.h"

// include implementations of utility functions
#include "geometry_utils.h"
#include "geometry_utils.inl"
#include "collision_utils.inl"
#include "exact_disk_utils.inl"
#include "diffusion_utils.inl"
#include "rxn_utils.inl"
#include "grid_utils.inl"
#include "wall_utils.inl"

using namespace std;
using namespace BNG;

namespace MCell {

void DiffuseReactEvent::step() {
  assert(world->get_partitions().size() == 1 && "Must extend cache to handle multiple partitions");
  assert(cmp_eq(event_time, (float_t)world->stats.get_current_iteration()) &&
  	"DiffuseReactEvent is expected to be run exactly at iteration starts");

  // for each partition
  for (Partition& p: world->get_partitions()) {
    // diffuse molecules that are scheduled for this iteration
    p.get_molecules_ready_for_diffusion(molecules_ready_array);
    diffuse_molecules(p, molecules_ready_array);
  }
}


void DiffuseReactEvent::diffuse_molecules(Partition& p, const MoleculeIdsVector& molecule_ids) {

  // we need to strictly follow the ordering in mcell3, therefore steps 2) and 3) do not use the time
  // for which they were scheduled but rather simply the order in which these "microevents" were created

  std::vector<DiffuseAction> delayed_diffusions;

#ifndef MCELL3_4_ALWAYS_SORT_MOLS_BY_TIME_AND_ID

  // 1) first diffuse already existing molecules
  uint existing_mols_count = molecule_ids.size();
  for (size_t i = 0; i < existing_mols_count; i++) {
    molecule_id_t id = molecule_ids[i];

    // here we compute both release delay and also partially sort molecules to be diffused
    // so that molecules not scheduled for he beginning of the event are diffused later
    // (but according to their id, not time)
    float_t diffusion_delay =  p.get_m(id).diffusion_time - event_time;

    if (cmp_eq(diffusion_delay, 0.0)) { // NOTE: this can be optimized
      // existing molecules or created at the beginning of this timestep
      // - simulate whole time step for this molecule
      diffuse_single_molecule(p, id, WallTileIndexPair());
    }
    else {
      // released during this iteration but not at the beginning, postpone its diffusion
      assert(diffusion_delay > 0 && diffusion_delay < time_up_to_next_barrier);
      delayed_diffusions.push_back(DiffuseAction(id));
    }
  }
#else
  for (molecule_id_t id: molecule_ids) {
    // merge with actions
    new_diffuse_actions.push_back(DiffuseAction(id));
  }
#endif

  // 2) mcell3 first handles diffusions of existing molecules, then the delayed diffusions
  // actions created by the diffusion of all these molecules are handled later
  for (size_t i = 0; i < delayed_diffusions.size(); i++) {
    const DiffuseAction& action = delayed_diffusions[i];

    Molecule& m = p.get_m(action.id);
    diffuse_single_molecule(
        p, action.id,
        action.where_created_this_iteration
    );
  }

  // 3) simulate remaining time of molecules created with reactions or
  // scheduled unimolecular reactions
  // need to call .size() each iteration because the size can increase,
  // again, we are using it as a queue and we do not follow the time when
  // they were created
#ifndef MCELL3_4_ALWAYS_SORT_MOLS_BY_TIME_AND_ID
  for (size_t i = 0; i < new_diffuse_actions.size(); i++) {

    DiffuseAction& action = new_diffuse_actions[i];
    diffuse_single_molecule(
        p, action.id,
        action.where_created_this_iteration // making a copy of this pair
    );
  }
#else
  for (size_t i = 0; i < new_diffuse_actions.size(); i++) {

    // find the closest event that we did not process yet
    uint best_index = UINT_INVALID;
    float_t best_time = TIME_FOREVER;
    for (size_t k = 0; k < new_diffuse_actions.size(); k++) {
      molecule_id_t id = new_diffuse_actions[k].id;
      if (id == MOLECULE_ID_INVALID) {
        continue;
      }
      const Molecule& m = p.get_m(id);
      if (m.diffusion_time < best_time) {
        best_index = k;
        best_time = m.diffusion_time;
      }
    }
    assert(best_index != UINT_INVALID);

    DiffuseAction& action = new_diffuse_actions[best_index];
    diffuse_single_molecule(
        p, action.id,
        action.where_created_this_iteration // making a copy of this pair
    );

    // invalidate this action
    new_diffuse_actions[best_index].id = MOLECULE_ID_INVALID;
  }
#endif

  new_diffuse_actions.clear();
}


inline float_t DiffuseReactEvent::get_max_time(Partition& p, const molecule_id_t m_id) {
  Molecule& m = p.get_m(m_id);
  const Species& species = p.get_species(m.species_id);

  float_t diffusion_time = m.diffusion_time;
  float_t unimol_rx_time = m.unimol_rx_time;
  float_t time_from_event_start = m.diffusion_time - event_time;

  // clamp to barrier time
  float_t max_time = time_up_to_next_barrier - time_from_event_start;

  // clamp to unimol_rx time
  if (unimol_rx_time != TIME_INVALID &&
      unimol_rx_time < diffusion_time + max_time) {
    assert(unimol_rx_time >= diffusion_time);
    max_time = unimol_rx_time - diffusion_time;
  }

  if (!m.has_flag(MOLECULE_FLAG_MATURE) &&
      species.time_step > DIFFUSE_REACT_EVENT_PERIODICITY) {
    // - newly created particles that have long time steps gradually increase
    //   their timestep to the full value,
    // - this behavior is here due to unbinding events so that the molecule
    //   does not diffuse so far when in reality it could re-bind
    float_t age_dependent_max_time = 1.0 + 0.2 * (diffusion_time - m.birthday);
    if (max_time > age_dependent_max_time) {
      max_time = age_dependent_max_time;
    }
    if (age_dependent_max_time > species.time_step) {
      // we are alive for long enough, no need to increase the time_step gradually anymore
      m.set_flag(MOLECULE_FLAG_MATURE);
    }
  }

  return max_time;
}


void DiffuseReactEvent::diffuse_single_molecule(
    Partition& p,
    const molecule_id_t m_id,
    WallTileIndexPair wall_tile_pair_where_created_this_iteration // set only for newly created molecules
) {
  Molecule& m_initial = p.get_m(m_id);
  float_t diffusion_start_time = m_initial.diffusion_time;
  assert(diffusion_start_time + EPS >= event_time && before_this_iterations_end(diffusion_start_time));
  assert(m_initial.birthday != TIME_INVALID && m_initial.birthday <= diffusion_start_time);

  if (m_initial.is_defunct()) {
    return;
  }

  if (m_initial.unimol_rx_time == diffusion_start_time) {
    // call to this diffuse_single_molecule was scheduled for time of an unimol rxn
    // may invalidate molecule references
    bool molecule_survived = react_unimol_single_molecule(p, m_id);
    if (!molecule_survived) {
      return;
    }
  }

  Molecule& m = p.get_m(m_id);
  // if the molecule is a "newbie", its unimolecular reaction was not yet scheduled,
  assert(
      !(m.has_flag(MOLECULE_FLAG_SCHEDULE_UNIMOL_RXN) && m.has_flag(MOLECULE_FLAG_RESCHEDULE_UNIMOL_RXN_ON_NEXT_RXN_RATE_UPDATE)) &&
      "Only one of these flags may be set"
  );
  if (m.has_flag(MOLECULE_FLAG_SCHEDULE_UNIMOL_RXN)) {
    m.clear_flag(MOLECULE_FLAG_SCHEDULE_UNIMOL_RXN);
    pick_unimol_rxn_class_and_set_rxn_time(p, diffusion_start_time, m);
  }

  // we might need to change the reaction rate right now
  if (m.has_flag(MOLECULE_FLAG_RESCHEDULE_UNIMOL_RXN_ON_NEXT_RXN_RATE_UPDATE) && cmp_eq(m.unimol_rx_time, diffusion_start_time)) {
    assert(m.unimol_rx_time != TIME_INVALID);

    m.clear_flag(MOLECULE_FLAG_RESCHEDULE_UNIMOL_RXN_ON_NEXT_RXN_RATE_UPDATE);
    pick_unimol_rxn_class_and_set_rxn_time(p, diffusion_start_time, m);
  }

  float_t unimol_rx_time = m.unimol_rx_time; // copy for a minor compiler optimization

#ifdef DEBUG_DIFFUSION
  const BNG::Species& debug_species = p.get_species(m.species_id);
  DUMP_CONDITION4(
#ifdef DUMP_NONDIFFUSING_VMS
    const char* title =
        (debug_species.can_diffuse()) ?
            (m.is_vol() ? "Diffusing vm:" : "Diffusing sm:") :
            (m.is_vol() ? "Not diffusing vm:" : "Not diffusing sm:");
    m.dump(p, "", title, world->get_current_iteration(), diffusion_start_time);
#else
    if (debug_species.can_diffuse()) {
      const char* title = (m.is_vol() ? "Diffusing vm:" : "Diffusing sm:");
      m.dump(p, "", title, world->get_current_iteration(), diffusion_start_time);
    }
#endif
  );
#endif

  // max_time is the time for which we should simulate the diffusion
  float_t max_time = get_max_time(p, m_id);
  
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

  // update time for which the molecule should be scheduled next
  Molecule& m_for_sched_update = p.get_m(m_id);
  if (!m_for_sched_update.is_defunct()) {
    float_t event_end_time = p.stats.get_current_iteration() + 1;

    const Species& species = world->get_all_species().get(m_for_sched_update.species_id);

    // TODO: for some reason, MCell3 calls diffusion of non-diffusible surface molecules,
    //       this is not needed since they cannot react but for compatibility we must keep this behavior
    if (species.can_diffuse() || species.has_flag(SPECIES_FLAG_CAN_SURFSURF)) {
      m_for_sched_update.diffusion_time += max_time;
      // - for MCell3 compatibility purposes, we must not keep any unimol rxn for the next iteration
      //   so we use precise comparison here
      // - on the other hand, we do not want to simulate diffusion for a tiny amount of time,
      //   so we use tolerance when checking whether we should keep diffusion itself for
      //   the next time, the error accumulation can be quite big, therefore we are using SQRT_EPS
      if (
          ( m_for_sched_update.unimol_rx_time != TIME_INVALID &&
              m_for_sched_update.unimol_rx_time < event_end_time
          ) ||
          before_this_iterations_end(m_for_sched_update.diffusion_time)
      ) {
        // reschedule molecule for this iteration because we did not use up all its time
        DiffuseAction diffuse_action(m_for_sched_update.id);
        new_diffuse_actions.push_back(diffuse_action);
      }
      else {
        // round diffusion_time to a whole number if it is close to it,
        // this did not give any error for the test that were available when this
        // code was implemented
        float_t rounded_dt = round_f(m_for_sched_update.diffusion_time);
        if (cmp_eq(m_for_sched_update.diffusion_time, rounded_dt, SQRT_EPS)) {
          m_for_sched_update.diffusion_time = rounded_dt;
        }
      }
    }
    else {
      // cannot diffuse
      if (m_for_sched_update.unimol_rx_time != TIME_INVALID) {
        // schedule for its unimol rxn
        m_for_sched_update.diffusion_time = m_for_sched_update.unimol_rx_time;

        // reschedule molecule for unimol rxn this iteration
        if (m_for_sched_update.unimol_rx_time < event_end_time) {
          DiffuseAction diffuse_action(m_for_sched_update.id);
          new_diffuse_actions.push_back(diffuse_action);
        }
      }
      else {
        // no need to schedule at all
        m_for_sched_update.diffusion_time = TIME_FOREVER;
      }
    }
  }
}

// ---------------------------------- volume diffusion ----------------------------------

void sort_collisions_by_time(CollisionsVector& molecule_collisions) {
  sort( molecule_collisions.begin(), molecule_collisions.end(),
      [ ]( const Collision& lhs, const Collision& rhs )
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
    float_t& max_time,
    const float_t diffusion_start_time,
    WallTileIndexPair& wall_tile_pair_where_created_this_iteration
) {
  Molecule& vm = p.get_m(vm_id);
  const BNG::Species& species = p.get_species(vm.species_id);

  if (!species.can_diffuse()) {
    return;
  }

  p.stats.inc_diffuse_3d_calls();

  // diffuse each molecule - get information on position change
  Vec3 displacement;

  float_t steps = 1.0;
  float_t t_steps = 1.0;
  float_t rate_factor = 1.0;
  float_t r_rate_factor = 1.0;
  DiffusionUtils::compute_vol_displacement(
      p, species, vm, max_time, world->rng,
      displacement, rate_factor, r_rate_factor, steps, t_steps
  );

#ifdef DEBUG_TIMING
  DUMP_CONDITION4(
      dump_vol_mol_timing(
          "- Timing vm", p.stats.get_current_iteration(), vm_id,
          diffusion_start_time, max_time, vm.unimol_rx_time,
          rate_factor, r_rate_factor, steps, t_steps
      );
  );
#endif

#ifdef DEBUG_DIFFUSION
  DUMP_CONDITION4(
      if (species.can_diffuse()) {
        displacement.dump("  displacement:", "");
        cout << "t_steps: " << t_steps << "\n";
      }
  );
#endif
  // note: we are ignoring use_expanded_list setting compared to mcell3

  // cut the displacement so that it does not cross partition?
  // or how to simplify this case when we know that we will definitely
  // be stopped by a wall (we need to fail if not anyway)

  Vec3 remaining_displacement = displacement;

  RayTraceState state;
  CollisionsVector molecule_collisions;
  bool was_defunct = false;
  wall_index_t last_hit_wall_index = WALL_INDEX_INVALID;

  float_t elapsed_molecule_time = diffusion_start_time; // == vm->t
  bool can_vol_react = species.can_vol_react();
  do {
    state =
        ray_trace_vol(
            p, world->rng,
            vm_id /* changes position */,
            can_vol_react,
            last_hit_wall_index,
            remaining_displacement,
            molecule_collisions
        );

    if (molecule_collisions.size() > 1) {
      sort_collisions_by_time(molecule_collisions);
    }

#ifdef DEBUG_COLLISIONS
    DUMP_CONDITION4(
        Collision::dump_array(p, molecule_collisions);
    );
#endif

    // evaluate and possible execute collisions and reactions
    for (Collision& collision: molecule_collisions) {

#if POS_T_BYTES == 8
      assert(collision.time >= 0 && collision.time <= 1);
#else
      assert(collision.time >= 0 && cmp_le(collision.time, 1, (float_t)POS_SQRT_EPS));
#endif

      if (collision.is_vol_mol_vol_mol_collision()) {
        p.stats.inc_vol_mol_vol_mol_collisions();

        // ignoring immediate collisions
        if (CollisionUtils::is_immediate_collision(collision.time)) {
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
        if (world->get_callbacks().needs_callback_for_mol_wall_hit(colliding_wall.object_id, vm_new_ref.species_id)) {
          // prepare part of information
          shared_ptr<API::MolWallHitInfo> info = make_shared<API::MolWallHitInfo>();
          info->molecule_id = vm_new_ref.id;
          info->geometry_object_id = colliding_wall.object_id; // I would need Model to be accessible here
          info->partition_wall_index = colliding_wall.index; // here as well
          info->time = elapsed_molecule_time + t_steps * collision.time;
          info->pos3d = collision.pos;
          info->time_before_hit = elapsed_molecule_time;
          info->pos3d_before_hit = vm_new_ref.v.pos;

          world->get_callbacks().do_mol_wall_hit_callback(info);
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
        if (p.get_species(vm_new_ref.species_id).has_flag(SPECIES_FLAG_CAN_VOLSURF) && colliding_wall.has_initialized_grid()) {
          int collide_res = collide_and_react_with_surf_mol(
              p, collision, r_rate_factor,
              wall_tile_pair_where_created_this_iteration,
              last_hit_wall_index,
              remaining_displacement,
              t_steps,
              elapsed_molecule_time
          );

          if (collide_res == 1) { // FIXME: use enum
            was_defunct = true;
            break;
          }
          else if (collide_res == 0) {
            // flip, position and counted vol change was already handled in outcome_products_random
            // continue with diffusion
            break;
          }
        }

        // check possible reaction with walls
        if (p.get_species(vm_new_ref.species_id).has_flag(SPECIES_FLAG_CAN_VOLWALL)) {
          WallRxnResult collide_res = collide_and_react_with_walls(p, collision, r_rate_factor, elapsed_molecule_time, t_steps);

          if (collide_res == WallRxnResult::Transparent) {
            // update molecules' counted volume, time and displacement and continue
            was_defunct = !cross_transparent_wall(
                p, collision, vm_new_ref, remaining_displacement,
                t_steps, elapsed_molecule_time, last_hit_wall_index
            );

            // continue with diffusion
            break;
          }
          else if (collide_res == WallRxnResult::Destroyed) {
            was_defunct = true;
          }
          else if (collide_res == WallRxnResult::Reflect) {
            // reflect in reflect_or_periodic_bc and continue with diffusion
          }
          else {
            assert(false);
          }
        }


        if (!was_defunct) {
          elapsed_molecule_time += t_steps * collision.time;
          // if a molecule was reflected, changes its position to the reflection point
          int res = CollisionUtils::reflect_or_periodic_bc(
              p, collision,
              vm_new_ref, remaining_displacement, t_steps, last_hit_wall_index
          );
          p.stats.inc_mol_wall_reflections();
          assert(res == 0 && "Periodic box BCs are not supported yet");
        }

        break; // we reflected and did not react, do ray_trace again
      }
    }

    assert(p.get_m(vm_id).v.subpart_index == p.get_subpart_index(p.get_m(vm_id).v.pos));

  } while (unlikely(state != RayTraceState::FINISHED && !was_defunct));

  if (!was_defunct) {
    // need to get a new reference
    Molecule& m_new_ref = p.get_m(vm_id);

    if (m_new_ref.is_vol()) {

      // change molecules' subpartition

#ifdef DEBUG_DIFFUSION
      DUMP_CONDITION4(
        // the subtraction of diffusion_time_step doesn't make much sense but is needed to make the dump the same as in mcell3
        // need to check it further
          m_new_ref.dump(p, "", "- Final vm position:", world->get_current_iteration(), 0);
      );
#endif

      // are we still in the same partition or do we need to move?
      bool move_to_another_partition = !p.in_this_partition(m_new_ref.v.pos);
      if (move_to_another_partition) {
        Vec3 pos_um = m_new_ref.v.pos * p.config.length_unit;
        const BNG::Species& s = p.get_species(m_new_ref.species_id);
        world->fatal_error(
            "Molecule with species " + s.name + " escaped the simulation area defined by partition size.\n"
            "Diffused molecule reached position (" +
            to_string(pos_um.x) + ", " + to_string(pos_um.y) + ", " + to_string(pos_um.z) + "). "
            "MCell4 requires a fixed-size simulation 3D space compared to MCell3 that allows unlimited space.\n"
            "One can create a geometry object box that keeps all molecules within a given area.\n"
            "This box can either reflect molecules (by default) or destroy the molecules with an absorptive surface class.\n"
            "Another option is to increase the partition size through CellBlender settings Partitions or\n"
            "through Model.config.partition_dimension."
        );
      }

      // change subpartition
      p.update_molecule_reactants_map(m_new_ref);
    }
  }
}


// collect possible collisions for molecule vm that has to displace by remaining_displacement,
// returns possible collisions in molecule_collisions, new position in new_pos and
// index of the new subparition in new_subpart_index
// later, this will check collisions until a wall is hit
// we assume that wall collisions do not occur so often
// inlining of this function does not help with performance
RayTraceState ray_trace_vol(
    Partition& p,
    rng_state& rng,
    const molecule_id_t vm_id, // molecule that we are diffusing, we are changing its pos  and possibly also subvolume
    const bool can_vol_react,
    const wall_index_t last_hit_wall_index, // is WALL_INDEX_INVALID when our molecule did not reflect from anything this diffusion step yet
    Vec3& remaining_displacement, // in/out - recomputed if there was a reflection
    CollisionsVector& collisions // both mol mol and wall collisions
    ) {
  Molecule& vm = p.get_m(vm_id);
  p.stats.inc_ray_voxel_tests();

  RayTraceState res_state = RayTraceState::FINISHED;
  collisions.clear();

  pos_t radius = p.config.rx_radius_3d;

  // if we would get out of this partition, cut off the displacement
  // so we check collisions just here
  Vec3 partition_displacement;
  if (!p.in_this_partition(vm.v.pos + remaining_displacement)) {
    partition_displacement = CollisionUtils::get_displacement_up_to_partition_boundary(p, vm.v.pos, remaining_displacement);
  }
  else {
    partition_displacement = remaining_displacement;
  }

  // first get what subpartitions might be relevant
  SubpartIndicesVector crossed_subparts_for_walls;
  SubpartIndicesSet crossed_subparts_for_molecules;

  CollisionUtils::collect_crossed_subparts(
      p, vm, partition_displacement,
      radius, p.config.subpartition_edge_length,
      can_vol_react, true,
      crossed_subparts_for_walls, crossed_subparts_for_molecules
  );

#ifndef NDEBUG
  if (can_vol_react) {
    // crossed subparts must contain our own subpart
    assert(crossed_subparts_for_molecules.count(vm.v.subpart_index) == 1);
    bool debug_found = false;
    for (subpart_index_t debug_index: crossed_subparts_for_walls) {
      if (debug_index == vm.v.subpart_index) {
        debug_found = true;
        break;
      }
    }
    assert(debug_found && "Did not find the starting subpartition in crossed subparts for walls");

    // also, each subpart from a wall must be present when checking molecules
    for (subpart_index_t debug_index: crossed_subparts_for_walls) {
      assert(crossed_subparts_for_molecules.count(debug_index) == 1);
    }
  }
#endif

  // changed when wall was hit
  Vec3 displacement_up_to_wall_collision = remaining_displacement;

  // check wall collisions in the crossed subparitions,
  if (!crossed_subparts_for_walls.empty()) {
    Collision closest_collision;
    for (subpart_index_t subpart_w_walls_index: crossed_subparts_for_walls) {

      bool collision_found =
          CollisionUtils::get_closest_wall_collision( // mcell3 does this only for the current subvol
              p,
              vm,
              subpart_w_walls_index,
              last_hit_wall_index,
              rng,
              remaining_displacement,
              displacement_up_to_wall_collision, // may be update in case we need to 'redo' the collision detection
              closest_collision
          );

      // stop at first crossing because crossed_subparts_for_walls are ordered
      // and we are sure that if we hit a wall in the actual subpartition, we cannot
      // possibly hit another wall in a subpartition that follows
      if (collision_found) {
        collisions.push_back(closest_collision);
        res_state = RayTraceState::RAY_TRACE_HIT_WALL;
        break;
      }
    }
  }

  if (can_vol_react) {
    if (res_state == RayTraceState::RAY_TRACE_HIT_WALL) {
      // recompute collect_crossed_subparts if there was a wall collision
      // NOTE: this can be in theory done more efficiently if we knew the order of subpartitions that we hit in the previous call
      crossed_subparts_for_molecules.clear();
      crossed_subparts_for_walls.clear();
      CollisionUtils::collect_crossed_subparts(
          p, vm, displacement_up_to_wall_collision,
          radius, p.config.subpartition_edge_length,
          true, false,
          crossed_subparts_for_walls, crossed_subparts_for_molecules
      );
    }

    // check molecule collisions for each SP
    for (subpart_index_t subpart_index: crossed_subparts_for_molecules) {
      // get cached reacting molecules for this SP
      const MoleculeIdsSet& sp_reactants = p.get_volume_molecule_reactants(subpart_index, vm.species_id);

      // for each molecule in this SP
      for (molecule_id_t colliding_vm_id: sp_reactants) {
        CollisionUtils::collide_mol_loop_body(
            p,
            vm,
            colliding_vm_id,
            remaining_displacement,// needs the full displacement to compute reaction time displacement_up_to_wall_collision,
            radius,
            collisions
        );
      }
    }
  }

  if (res_state == RayTraceState::FINISHED) {
    vm.v.pos = vm.v.pos + remaining_displacement;
    vm.v.subpart_index = p.get_subpart_index(vm.v.pos);
  }

  return res_state; // no wall was hit
}


// handle collision of two volume molecules: checks probability of reaction,
// executes this reaction, removes reactants and creates products
// returns true if the first reactant was destroyed
bool DiffuseReactEvent::collide_and_react_with_vol_mol(
    Partition& p,
    const Collision& collision,
    Vec3& displacement,
    const float_t t_steps,
    const float_t r_rate_factor,
    const float_t elapsed_molecule_time
)  {

  Molecule& colliding_molecule = p.get_m(collision.colliding_molecule_id); // am
  Molecule& diffused_molecule = p.get_m(collision.diffused_molecule_id); // m

  // returns 1 when there are no walls at all
  float_t factor = ExactDiskUtils::exact_disk(
      p, collision.pos, displacement, p.config.rx_radius_3d,
      diffused_molecule, colliding_molecule,
      p.config.use_expanded_list
  );

  if (factor < 0) {
    return 0; // reaction blocked by a wall
  }

  RxnClass* rxn_class = collision.rxn_class;
  assert(rxn_class != nullptr);

  float_t absolute_collision_time = elapsed_molecule_time + t_steps * collision.time;

  //  rx->prob_t is always NULL in out case update_probs(world, rx, m->t);
  // returns which reaction pathway to take
  float_t scaling = factor * r_rate_factor;
  int i = RxnUtils::test_bimolecular(
      rxn_class, world->rng, colliding_molecule, diffused_molecule, scaling, 0, absolute_collision_time);

  if (i < RX_LEAST_VALID_PATHWAY) {
    return false;
  }
  else {
    // might invalidate references
    int j = outcome_bimolecular(p, collision, i, absolute_collision_time);
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
// TODO: return enum
int DiffuseReactEvent::collide_and_react_with_surf_mol(
    Partition& p,
    const Collision& collision,
    const float_t r_rate_factor,
    WallTileIndexPair& where_created_this_iteration,
    wall_index_t& last_hit_wall_index,
    Vec3& remaining_displacement,
    float_t& t_steps,
    float_t& elapsed_molecule_time
) {
  Wall& wall = p.get_wall(collision.colliding_wall_index);
  Grid& grid = wall.grid;

  tile_index_t j = GridUtils::xyz2grid_tile_index(p, collision.pos, wall);
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
  RxnUtils::trigger_bimolecular(
    p.bng_engine,
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

  float_t collision_time = elapsed_molecule_time + t_steps * collision.time;

  int selected_rx_pathway;
  if (matching_rxn_classes.size() == 1) {
    selected_rx_pathway = RxnUtils::test_bimolecular(
        matching_rxn_classes[0], world->rng,
        diffused_molecule, colliding_molecule,
        scaling_coefs[0], 0, collision_time);
  }
  else {
    mcell_error("Internal error: multiple rxn classes should not be needed for surf mol rxn.");
  }

  if (selected_rx_pathway < RX_LEAST_VALID_PATHWAY) {
    return -1; /* No reaction */
  }

  /* run the reaction */
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
  else if (outcome_bimol_result == RX_FLIP) {
    // update position and timing for flipped molecule

    float_t t_smash = collision.time;
    Molecule& vm = p.get_m(collision.diffused_molecule_id);

    // update position and subpart if needed
    vm.v.pos = collision.pos;
    vm.v.subpart_index = p.get_subpart_index(vm.v.pos);

    CollisionUtils::update_counted_volume_id_when_crossing_wall(
        p, wall, collision.get_orientation_against_wall(), vm);

    // TODO: same code is on multiple places, e.g. in cross_transparent_wall,
    // make a function for it
    remaining_displacement = remaining_displacement * Vec3(1.0 - t_smash);
    elapsed_molecule_time += t_steps * t_smash;

    t_steps *= (1.0 - t_smash);
    if (t_steps < EPS) {
      t_steps = EPS;
    }

    last_hit_wall_index = wall.index;
    return 0;
  }
  else {
    last_hit_wall_index = wall.index;
    return -1;
  }
}


/******************************************************************************
 *
 * check_collisions_with_walls is a helper function used in diffuse_3D to handle
 * collision of a diffusing molecule with a wall
 *
 * Return values:
 *
 * -1 : nothing happened - continue on with next smash targets
 *  0 : reaction happened and we still exist but are done with the current smash
 *      target
 *  1 : reaction happened and we are destroyed
 *
 ******************************************************************************/
inline WallRxnResult DiffuseReactEvent::collide_and_react_with_walls(
    Partition& p,
    Collision& collision,
    const float_t r_rate_factor,
    const float_t elapsed_molecule_time,
    const float_t t_steps
) {
  assert(collision.type == CollisionType::WALL_FRONT || collision.type == CollisionType::WALL_BACK);

  Molecule& diffused_molecule = p.get_m(collision.diffused_molecule_id); // m
  assert(diffused_molecule.is_vol());

  Wall& wall = p.get_wall(collision.colliding_wall_index);

  RxnClassesVector matching_rxn_classes;
  orientation_t orient = (collision.type == CollisionType::WALL_FRONT) ? ORIENTATION_UP : ORIENTATION_DOWN;
  RxnUtils::trigger_intersect(p, diffused_molecule, orient, wall, true, matching_rxn_classes);
  if (matching_rxn_classes.empty() ||
      (matching_rxn_classes.size() == 1 && matching_rxn_classes[0]->type == BNG::RxnType::Reflect)) {
    return WallRxnResult::Reflect;
  }

  const BNG::RxnClass* transp_rxn_class = nullptr;
  for (const BNG::RxnClass* rxn: matching_rxn_classes) {
    if (rxn->is_transparent()) {
      transp_rxn_class = rxn;
      break;
    }
  }

  if (transp_rxn_class != nullptr) {
    // we are crossing this wall
    assert(matching_rxn_classes.size() == 1 && "Expected only a single transparent rxn for transparent walls");
    return WallRxnResult::Transparent;
  }

  /* Collisions with the surfaces declared REFLECTIVE are treated similar to
   * the default surfaces after this loop. */

  // TODO: use rxn_class_index_t and rxn_index_t also in other places
  rxn_class_index_t rxn_class_index = RNX_CLASS_INDEX_INVALID;
  rxn_class_pathway_index_t pathway_index;

  float_t current_time = elapsed_molecule_time + t_steps * collision.time;

  if (matching_rxn_classes.size() == 1) {
    rxn_class_index = 0;
    pathway_index = RxnUtils::test_intersect(matching_rxn_classes[0], r_rate_factor, current_time, world->rng);
  } else {
    pathway_index = RxnUtils::test_many_intersect(
        matching_rxn_classes, r_rate_factor, current_time, rxn_class_index, world->rng);
  }

  if (rxn_class_index != RNX_CLASS_INDEX_INVALID && pathway_index >= PATHWAY_INDEX_LEAST_VALID) {
    BNG::RxnClass* rxn_class = matching_rxn_classes[rxn_class_index];

    assert(collision.type == CollisionType::WALL_FRONT || collision.type == CollisionType::WALL_BACK);

    int j = outcome_intersect(
        p, rxn_class, pathway_index, collision, current_time);

    if (j == RX_FLIP) {
      return WallRxnResult::Transparent;
    }
    else if (j == RX_A_OK) {
      return WallRxnResult::Reflect;
    }
    else if (j == RX_DESTROY) {
      return WallRxnResult::Destroyed;
    }
    else {
      assert(false);
      return WallRxnResult::Invalid;
    }
  }

  return WallRxnResult::Reflect;
}


// ---------------------------------- surface diffusion ----------------------------------
inline void DiffuseReactEvent::diffuse_surf_molecule(
    Partition& p,
    const molecule_id_t sm_id,
    float_t& max_time,
    const float_t diffusion_start_time
) {
  Molecule& sm = p.get_m(sm_id);
  const BNG::Species& species = p.get_species(sm.species_id);

  float_t steps = 0.0;
  float_t t_steps = 0.0;

  wall_index_t original_wall_index = sm.s.wall_index;

  /* Where are we going? */
  if (species.get_time_step() > max_time) {
    t_steps = max_time;
  }
  else {
    t_steps = species.get_time_step();
  }

  bool diffusible = species.can_diffuse();

  if (diffusible) {
    if (species.get_time_step() > max_time) {
      steps = max_time / species.get_time_step() ;
      if (steps < EPS) {
        t_steps = EPS * species.get_time_step();
        steps = EPS;
      }
    }
    else {
      steps = 1.0;
    }

    float_t space_factor = 0.0;

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

      Vec2 displacement;
      DiffusionUtils::compute_surf_displacement(species, space_factor, world->rng, displacement);

#ifdef DEBUG_DIFFUSION
      DUMP_CONDITION4(
          if (species.can_diffuse()) {
            displacement.dump("  displacement:", "");
          }
      );
#endif

      assert(!species.has_flag(SPECIES_FLAG_SET_MAX_STEP_LENGTH) && "not supported yet");

      // ray_trace does the movement and all other stuff
      Vec2 new_loc;
      BNG::RxnClass* absorb_now_rxn_class;
      wall_index_t new_wall_index =
          ray_trace_surf(p, species, sm_id, displacement, new_loc, absorb_now_rxn_class);

      // Either something ambiguous happened or we hit absorptive border
      if (new_wall_index == WALL_INDEX_INVALID) {
        if (absorb_now_rxn_class != nullptr) {
          absorb_now_rxn_class->init_rxn_pathways_and_rates();
          assert(absorb_now_rxn_class->get_num_pathways() == 1);
          bool kept = outcome_unimolecular(
              p, p.get_m(sm_id), diffusion_start_time, absorb_now_rxn_class, 0, nullptr);
          assert(!kept);
          return;
        }

        continue; /* Something went wrong--try again */
      }

      assert(absorb_now_rxn_class == nullptr);

      // After diffusing, are we still on the SAME triangle?
      if (new_wall_index == sm.s.wall_index) {
        if (!DiffusionUtils::move_sm_on_same_triangle(p, sm, new_loc)) {
          // we must try a different position
          continue;
        }
      }
      // After diffusing, we ended up on a NEW triangle.
      else {
        if (!DiffusionUtils::move_sm_to_new_triangle(p, sm, new_loc, new_wall_index)) {
          // we must try a different position
          continue;
        }

        // Check if we shouldn't reschedule the unimol rxn defined by
        // a surf mol + surf class rxn
        // Do this only if the current unimol rxn is not scheduled
        // for this timestep (t_steps is set to end at the unimol rx time)
        float_t time_until_unimol = sm.unimol_rx_time - t_steps - sm.diffusion_time;
        time_until_unimol = (time_until_unimol < 0) ? 0 : time_until_unimol;

        if (species.has_flag(SPECIES_FLAG_CAN_SURFWALL) &&
            (sm.unimol_rx_time == TIME_INVALID ||
             // using this overly complex condition because of MCell3 compatibility
             (time_until_unimol > EPS || time_until_unimol > EPS * (sm.diffusion_time + t_steps))
            )
        ) {
          sm.unimol_rx_time = TIME_INVALID;
          sm.set_flag(MOLECULE_FLAG_SCHEDULE_UNIMOL_RXN);
        }
      }
      find_new_position = 0;
    }
  } // if (species.space_step != 0)


  // TODO: what about molecules that cannot diffuse? they should not react,
  // but in MCell3 they do
  bool sm_still_exists = true;
  assert(!species.has_flag(SPECIES_FLAG_CAN_SURFSURFSURF) && "Not supported");
  bool can_surf_surf_react = species.has_flag(SPECIES_FLAG_CAN_SURFSURF);
  if (can_surf_surf_react && !species.is_target_only()) {
    assert(!species.has_flag(SPECIES_MOL_FLAG_TARGET_ONLY) && "Not sure what to do here");
    // the time t_steps should tell when the reaction occurred and it is quite weird because
    // it has nothing to do with the time spent diffusing
    sm_still_exists = react_2D_all_neighbors(p, sm, t_steps, diffusion_start_time);
  }

  // NOTE: surf-surf reactions on the same wall have higher priority,
  // this must be randomized once intermembrane rxns will be more used
  if (sm_still_exists && species.has_flag(SPECIES_FLAG_CAN_INTERMEMBRANE_SURFSURF) && !species.is_target_only()) {
    sm_still_exists = react_2D_intermembrane(p, sm, t_steps, diffusion_start_time);
  }

  if (sm_still_exists) {
    // reactions in react_2D_all_neighbors could have invalidated the molecules array
    Molecule& new_m_ref = p.get_m(sm_id);

    // for some reason, mcell3 defines a new unimol time if the molecule has moved
    bool changed_wall = new_m_ref.s.wall_index != original_wall_index;

    // we don't have to remove the molecule from the schedule, we can just change its unimol_rx_time,
    // this time is checked and against the scheduled time
    // mcell3 compatibility: we might change the schedule only if it is not already scheduled for this time step
    // NOTE: this condition looks a bit weird, some explanation would be useful
    if ((diffusible || can_surf_surf_react) &&
        ((!diffusible || changed_wall) &&
          new_m_ref.unimol_rx_time >= event_time + time_up_to_next_barrier)) { // ?? is usage or barrier_time_from_event_time correct?
      new_m_ref.unimol_rx_time = TIME_INVALID;
      new_m_ref.set_flag(MOLECULE_FLAG_SCHEDULE_UNIMOL_RXN);
    }
  }

  // update max_time - for how long the molecule diffused
  max_time = t_steps;

  // and log stats
  p.stats.inc_diffusion_cummtime(steps);
}


// returns true if molecule survived
bool DiffuseReactEvent::react_2D_all_neighbors(
    Partition& p,
    Molecule& sm,
    const float_t time, // same argument as t passed in mcell3 (come up with a better name)
    const float_t diffusion_start_time // diffusion_start_time + elapsed_molecule_time should be the time when reaction occurred
) {

#ifdef DEBUG_TIMING
  DUMP_CONDITION4(
      dump_react_2D_all_neighbors_timing(
          time, diffusion_start_time + elapsed_molecule_time
      );
  );
#endif

  const Wall& wall = p.get_wall(sm.s.wall_index);

  TileNeighborVector neighbors;
  GridUtils::find_neighbor_tiles(p, &sm, wall, sm.s.grid_tile_index, false, true, neighbors);

  if (neighbors.empty()) {
    return true;
  }

  const BNG::Species& sm_species = p.get_species(sm.species_id);

  // each item in array corresponds to one potential reaction
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
    const BNG::Species& nsm_species = p.get_species(nsm.species_id);

#ifdef DEBUG_RXNS
    DUMP_CONDITION4(
      // the subtraction of diffusion_time_step doesn't make much sense but is needed to make the dump the same as in mcell3
      // need to check it further
      nsm.dump(p, "", "  checking in react_2D_all_neighbors: ", world->get_current_iteration(), 0.0);
    );
#endif

    /* check whether the neighbor molecule is behind
       the restrictive region boundary   */
    if ((sm_species.has_flag(SPECIES_FLAG_CAN_REGION_BORDER) || nsm_species.has_flag(SPECIES_FLAG_CAN_REGION_BORDER)) &&
        sm.s.wall_index != nsm.s.wall_index
    ) {
      /* INSIDE-OUT check */
      if (WallUtils::walls_belong_to_at_least_one_different_restricted_region(p, wall, sm, nwall, nsm)) {
        continue;
      }

      /* OUTSIDE-IN check */
      // note: the pairing wall is same as in mcell3, TODO: explain why is it so
      if (WallUtils::walls_belong_to_at_least_one_different_restricted_region(p, wall, nsm, nwall, sm)) {
        continue;
      }
    }

    // returns value >=1 if there can be a reaction
    size_t orig_num_rxsn = matching_rxn_classes.size();
    RxnUtils::trigger_bimolecular_orientation_from_mols(
        p.bng_engine,
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

  size_t num_matching_rxn_classes = matching_rxn_classes.size();
  if (num_matching_rxn_classes == 0) {
    return true;
  }

  rxn_class_pathway_index_t selected_pathway_index;
  Collision collision;

  float_t collision_time = diffusion_start_time;

  /* Calculate local_prob_factor for the reaction probability.
     Here we convert from 3 neighbor tiles (upper probability
     limit) to the real "num_nbrs" neighbor tiles. */
  float_t local_prob_factor = 3.0 / neighbors.size();
  int rxn_class_index;
  if (num_matching_rxn_classes == 1) {
    // figure out what should happen
    selected_pathway_index = RxnUtils::test_bimolecular(
        matching_rxn_classes[0], world->rng,
        sm, p.get_m(reactant_molecule_ids[0]),
        correction_factors[0], local_prob_factor, collision_time);

    // there is just one possible class == one pair of reactants
    rxn_class_index = 0;
  }
  else {
    bool all_neighbors_flag = true;
    rxn_class_index =
        RxnUtils::test_many_bimolecular(
            matching_rxn_classes, correction_factors, local_prob_factor,
            world->rng, all_neighbors_flag, collision_time, selected_pathway_index);
    selected_pathway_index = 0; // TODO_PATHWAYS: use value from test_many_bimolecular
  }

  if (rxn_class_index == PATHWAY_INDEX_NO_RXN || selected_pathway_index < PATHWAY_INDEX_LEAST_VALID) {
    return true; /* No reaction */
  }

  collision = Collision(
      CollisionType::SURFMOL_SURFMOL,
      &p, sm.id, collision_time,
      reactant_molecule_ids[rxn_class_index],
      matching_rxn_classes[rxn_class_index]
  );

  /* run the reaction */
  int outcome_bimol_result = outcome_bimolecular(
      p, collision, selected_pathway_index, collision_time
  );

  return outcome_bimol_result != RX_DESTROY;
}


bool DiffuseReactEvent::react_2D_intermembrane(
    Partition& p, Molecule& sm,
    const float_t t_steps, const float_t diffusion_start_time
) {
  const Wall& w1 = p.get_wall(sm.s.wall_index);
  const Vec3& w1_vert0 = p.get_wall_vertex(w1, 0);

  // get 3d position
  Vec3 reac1_pos3d = GeometryUtils::uv2xyz(sm.s.pos, w1, w1_vert0);

  // which subpart are we in?
  // we check walls whose part is in the current subpartition, not the neighbors
  IVec3 subpart_indices;
  p.get_subpart_3d_indices(reac1_pos3d, subpart_indices);
  subpart_index_t subpart_index = p.get_subpart_index_from_3d_indices(subpart_indices);

  // with what species we may react? (this should be cached)
  const BNG::SpeciesRxnClassesMap* rxns_classes_map =
      p.get_all_rxns().get_bimol_rxns_for_reactant(sm.species_id);
  if (rxns_classes_map == nullptr || rxns_classes_map->empty()) {
    assert(false && "No rxns");
    return true;
  }

  uint_set<species_id_t> reacting_species;
  for (const auto& second_reactant_info: *rxns_classes_map) {
    const BNG::RxnClass* rxn_class = second_reactant_info.second;

    if (!rxn_class->is_intermembrane_surf_surf_rxn_class()) {
      continue;
    }

    species_id_t second_species_id = second_reactant_info.first;
    reacting_species.insert(second_species_id);
  }

  pos_t rxn_radius2 = p.config.intermembrane_rx_radius_3d * p.config.intermembrane_rx_radius_3d;

  typedef pair<molecule_id_t, pos_t> IdDist2Pair;
  small_vector<IdDist2Pair> reactants_dist2;

  // get neighboring subparts - this is necessary because
  // subpartitioning can put a boundary right between membranes
  SubpartIndicesSet subpart_indices_set;
  CollisionUtils::collect_neighboring_subparts(
      p, reac1_pos3d, subpart_indices, p.config.intermembrane_rx_radius_3d, p.config.subpartition_edge_length,
      subpart_indices_set
  );
  // and include current subpart
  subpart_indices_set.insert(subpart_index);

  // collect walls (this can be also somehow pre-computed)
  uint_set<wall_index_t> wall_indices;
  for (subpart_index_t si: subpart_indices_set) {
    // for each wall in this subpart that belongs to a different object
    const WallsInSubpart& subpart_wall_indices = p.get_subpart_wall_indices(subpart_index);
    wall_indices.insert(subpart_wall_indices.begin(), subpart_wall_indices.end());
  }

  for (wall_index_t wi: wall_indices) {
    const Wall& w2 = p.get_wall(wi);
    if (w2.object_id == w1.object_id) {
      continue;
    }

    const Grid& g = w2.grid;
    if (!g.is_initialized() || g.get_num_occupied_tiles() == 0) {
      continue;
    }
    const Vec3& w2_vert0 = p.get_wall_vertex(w2, 0);

    // not really efficient, goes through all tiles
    small_vector<molecule_id_t> molecule_ids;
    g.get_contained_molecules(molecule_ids);
    assert(molecule_ids.size() == g.get_num_occupied_tiles());

    for (molecule_id_t reac2_id: molecule_ids) {
      const Molecule& reac2 = p.get_m(reac2_id);
      assert(reac2.is_surf());
      if (reacting_species.count(reac2.species_id) == 0) {
        continue;
      }

      // compute distance
      Vec3 reac2_pos3d = GeometryUtils::uv2xyz(reac2.s.pos, w2, w2_vert0);
      pos_t dist2 = len3_squared(reac1_pos3d - reac2_pos3d);
      if (dist2 > rxn_radius2) {
        continue;
      }

      reactants_dist2.push_back(make_pair(reac2_id, dist2));
    }
  }

  if (reactants_dist2.empty()) {
    return true;
  }

  // sort potential reactants by distance
  sort(reactants_dist2.begin(), reactants_dist2.end(),
      [ ]( const IdDist2Pair& lhs, const IdDist2Pair& rhs ) {
        return lhs.second < rhs.second;
      }
  );

  // and similarly as in diffuse_volume_molecule, evaluate each potential collision

  // collision_time is the same as for surf mols,
  // TODO: not sure if correct, this molecule has already diffused and newly created
  // molecules birthdays will be this one
  float_t collision_time = diffusion_start_time;

  bool was_defunct = false;
  for (const IdDist2Pair& id_dist2_pair: reactants_dist2) {
    Molecule& sm2 = p.get_m(id_dist2_pair.first);

    // get what reaction should happen, the default orientation of molecules is UP and
    // also the rxn rule's pattern expects UP
    RxnClassesVector matching_rxn_classes;
    RxnUtils::trigger_bimolecular(
        p.bng_engine, sm, sm2, sm.s.orientation, sm2.s.orientation, matching_rxn_classes);

    if (matching_rxn_classes.empty()) {
      assert(false && "We already filtered-out molecules that can react");
      continue;
    }

    int selected_pathway_index = RX_NO_RX;
    int rxn_class_index;
    if (matching_rxn_classes.size() == 1) {
      selected_pathway_index = RxnUtils::test_bimolecular(
          matching_rxn_classes[0], world->rng,
          sm, sm2,
          // TODO: not sure what to put as scaling coeff
          1,
          0, // local prob factor
          collision_time);

        // there is just one possible class == one pair of reactants
        rxn_class_index = 0;
    }
    else {
      release_assert(false && "Multiple rxn classes for intermembrane rxns are not supported yet");
    }

    if (selected_pathway_index < RX_LEAST_VALID_PATHWAY) {
      continue; // No reaction
    }

    Collision collision = Collision(
        CollisionType::INTERMEMBRANE_SURFMOL_SURFMOL,
        &p, sm.id, collision_time,
        sm2.id,
        matching_rxn_classes[rxn_class_index]
    );

    /* run the reaction */
    int outcome_bimol_result = outcome_bimolecular(
        p, collision, selected_pathway_index, collision_time
    );
  }
  return !was_defunct;
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
    const BNG::Species& species,
    const molecule_id_t sm_id,
    Vec2& remaining_displacement,
    Vec2& new_pos,
    BNG::RxnClass*& absorb_now_rxn_class // set to non-null is molecule has to be absorbed
) {
  absorb_now_rxn_class = nullptr;

  Molecule& sm  = p.get_m(sm_id);
  const Wall* this_wall = &p.get_wall(sm.s.wall_index);

  Vec2 orig_pos = sm.s.pos;
  Vec2 this_pos = sm.s.pos;
  Vec2 this_disp = remaining_displacement;

  /* Will break out with return or break when we're done traversing walls */
  while (1) {

    int this_wall_edge_region_border = 0;
    //bool absorb_now = 0;

    /* Index of the wall edge that the SM hits */
    Vec2 boundary_pos;
    edge_index_t edge_index_that_was_hit =
        GeometryUtils::find_edge_point(*this_wall, this_pos, this_disp, boundary_pos);

    // Ambiguous edge collision. Give up and try again from diffuse_2D.
    if (edge_index_that_was_hit == EDGE_INDEX_CANNOT_TELL) {
      sm.s.pos = orig_pos;
      return WALL_INDEX_INVALID;
    }

    // We didn't hit the edge. Stay inside this wall. We're done!
    else if (edge_index_that_was_hit == EDGE_INDEX_WITHIN_WALL) {
      new_pos = this_pos + this_disp;
      sm.s.pos = orig_pos;
      return this_wall->index;
    }

    // Neither ambiguous (EDGE_INDEX_CANNOT_TELL) nor inside wall (EDGE_INDEX_WITHIN_WALL),
    // must have hit edge (0, 1, 2)
    Vec2 old_pos = this_pos;

    /* We hit the edge - check for the reflection/absorption from the
       edges of the wall if they are region borders
       Note - here we test for potential collisions with the region
       border while moving INSIDE OUT */
    bool reflect_now = false;
    if (species.can_interact_with_border()) {
      DiffusionUtils::reflect_absorb_inside_out(
          p, sm, *this_wall, edge_index_that_was_hit,
          reflect_now, absorb_now_rxn_class
      );

      if (absorb_now_rxn_class != nullptr) {
        return WALL_INDEX_INVALID;
      }
    }

    /* no reflection - keep going */
    if (!reflect_now) {
      wall_index_t target_wall_index =
          GeometryUtils::traverse_surface(*this_wall, old_pos, edge_index_that_was_hit, this_pos);

      if (target_wall_index != WALL_INDEX_INVALID) {
        /* We hit the edge - check for the reflection/absorption from the
           edges of the wall if they are region borders
           Note - here we test for potential collisions with the region
           border while moving OUTSIDE IN */
        if (species.can_interact_with_border()) {

          const Wall& target_wall = p.get_wall(target_wall_index);
          DiffusionUtils::reflect_absorb_outside_in(
              p, sm, target_wall, *this_wall,
              reflect_now, absorb_now_rxn_class
          );

          if (absorb_now_rxn_class != nullptr) {
            return WALL_INDEX_INVALID;
          }
        }

        if (!reflect_now) {
          this_disp = old_pos + this_disp;

#ifndef NDEBUG
          Edge& e = const_cast<Edge&>(this_wall->edges[edge_index_that_was_hit]);
          assert(e.is_initialized());
          e.debug_check_values_are_uptodate(p);
#endif

          Vec2 tmp_disp;
          GeometryUtils::traverse_surface(*this_wall, this_disp, edge_index_that_was_hit, tmp_disp);
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
    Vec2 new_disp = this_disp - (boundary_pos - old_pos);

    switch (edge_index_that_was_hit) {
      case EDGE_INDEX_0:
        new_disp.v *= -1.0;
        break;
      case EDGE_INDEX_1: {
        pos_t f;
        Vec2 reflector;
        reflector.u = -this_wall->uv_vert2.v;
        reflector.v = this_wall->uv_vert2.u - this_wall->uv_vert1_u;
        f = 1.0 / sqrt_p( len2_squared(reflector) );
        reflector *= f;
        f = 2.0 * dot2(new_disp, reflector);
        new_disp -= Vec2(f) * reflector;
        break;
      }
      case EDGE_INDEX_2: {
        pos_t f;
        Vec2 reflector;
        reflector.u = this_wall->uv_vert2.v;
        reflector.v = -this_wall->uv_vert2.u;
        f = 1.0 / sqrt_p( len2_squared(reflector) );
        reflector *= f;
        f = 2.0 * dot2(new_disp, reflector);
        new_disp -= Vec2(f) * reflector;
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


void DiffuseReactEvent::pick_unimol_rxn_class_and_set_rxn_time(
    Partition& p,
    const float_t current_time,
    Molecule& m
) {
  assert(current_time >= 0);

  BNG::RxnClassesVector rxn_classes;
  RxnUtils::pick_unimol_rxn_classes(p, m, current_time, rxn_classes);
  if (rxn_classes.empty()) {
    m.unimol_rx_time = TIME_INVALID;
    return;
  }

  uint idx = 0;
  if (rxn_classes.size() > 1) {
    idx = RxnUtils::test_many_unimol(rxn_classes, world->rng);
  }

  // there is a check when computing the reaction rate to make sure that the time is reasonably higher than 0
  float_t time_from_now = RxnUtils::compute_unimol_lifetime(p, world->rng, rxn_classes[idx], current_time, m);

  float_t scheduled_time = current_time + time_from_now;

  // we need to store the end time to the molecule because it is needed in diffusion to
  // figure out whether we should do the whole time step
  m.unimol_rx_time = scheduled_time;
}


// based on mcell3's check_for_unimolecular_reaction
// might invalidate vm references
// returns true if molecule survived and should be diffused right away
bool DiffuseReactEvent::react_unimol_single_molecule(
    Partition& p,
    const molecule_id_t m_id
) {
  // the unimolecular reaction class was already selected
  Molecule& m = p.get_m(m_id);

  if (m.is_defunct()) {
    return false;
  }

  float_t scheduled_time = m.unimol_rx_time;

  assert(!m.has_flag(MOLECULE_FLAG_SCHEDULE_UNIMOL_RXN));
  assert(scheduled_time >= event_time && scheduled_time <= event_time + time_up_to_next_barrier);

  if (m.has_flag(MOLECULE_FLAG_RESCHEDULE_UNIMOL_RXN_ON_NEXT_RXN_RATE_UPDATE)) {
    // this time event does not mean to execute the unimol action, but instead to
    // update the scheduled time for the updated reaction rate
    m.clear_flag(MOLECULE_FLAG_RESCHEDULE_UNIMOL_RXN_ON_NEXT_RXN_RATE_UPDATE);

    pick_unimol_rxn_class_and_set_rxn_time(p, scheduled_time, m);

    return true;
  }
  else {
    BNG::RxnClassesVector rxn_classes;
    RxnUtils::pick_unimol_rxn_classes(p, m, scheduled_time, rxn_classes);
    for (RxnClass* rxn_class: rxn_classes) {
      rxn_class->update_rxn_rates_if_needed(scheduled_time);
    }

    if (rxn_classes.empty()) {
      // there is no unimol rxn, this can happen for instance when
      // a surface molecule moved out of a regions with a surface class,
      // we must recompute the unimol lifetime
      pick_unimol_rxn_class_and_set_rxn_time(p, scheduled_time, m);
      return true;
    }

    uint idx = 0;
    if (rxn_classes.size() > 1) {
      idx = RxnUtils::test_many_unimol(rxn_classes, world->rng);
    }

    rxn_class_pathway_index_t pi = RxnUtils::which_unimolecular(m, rxn_classes[idx], world->rng);

    if (rxn_classes[idx]->is_unimol()) {
      // standard unimolecular rxn
      return outcome_unimolecular(p, m, scheduled_time, rxn_classes[idx], pi);
    }
    else {
      // surf mol + surf class rxn
      release_assert(m.is_surf() && "Only for surf mol + surf class rxns");
      const Wall& w = p.get_wall(m.s.wall_index);
      const Vec3& v0 = p.get_wall_vertex(w, 0);
      // orientation front or back is not important
      Vec3 pos = GeometryUtils::uv2xyz(m.s.pos, w, v0);
      Collision collision(CollisionType::WALL_FRONT, &p, m.id, scheduled_time, pos, w.index);
      return outcome_intersect(p, rxn_classes[idx], pi, collision, scheduled_time);
    }
  }
}


// checks if reaction should probabilistically occur and if so,
// destroys reactants
// returns RX_DESTROY when the primary reactant was destroyed, RX_A_OK if the reactant A was kept
// TODO: change return type for enum
int DiffuseReactEvent::outcome_bimolecular(
    Partition& p,
    const Collision& collision,
    const int path,
    const float_t time
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

  if (result == RX_A_OK || result == RX_FLIP) {
    Molecule& reacA = p.get_m(collision.diffused_molecule_id);
    Molecule& reacB = p.get_m(collision.colliding_molecule_id);

#ifdef DEBUG_RXNS
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
      return result;
    }
    else {
      return RX_DESTROY;
    }
  }

  return result;
}


/*************************************************************************
outcome_intersect:
  In: world: simulation state
      rx: reaction that's taking place
      path: path the reaction's taking
      surface: wall that is being struck
      reac: molecule that is hitting the wall
      orient: orientation of the molecule
      t: time that the reaction is occurring
      hitpt: location of collision with wall
      loc_okay:
  Out: Value depending on outcome:
       RX_A_OK if the molecule reflects
       RX_FLIP if the molecule passes through
       RX_DESTROY if the molecule stops, is destroyed, etc.
       Additionally, products are created as needed.
  Note: Can assume molecule is always first in the reaction.
*************************************************************************/
int DiffuseReactEvent::outcome_intersect(
    Partition& p,
    RxnClass* rxn_class,
    const rxn_class_pathway_index_t pathway_index,
    Collision& collision, // rxn_class can be set
    const float_t time
) {
  if (!rxn_class->is_standard()) {
    if (rxn_class->is_reflect()) {
      return RX_A_OK; /* just reflect */
    }
    else {
      assert(rxn_class->is_transparent()); // already dealt with before for better performance, but let's keep this general
      return RX_FLIP; /* Flip = transparent is default special case */
    }
  }

  /* If reaction object has ALL_MOLECULES or ALL_VOLUME_MOLECULES as the
   * first reactant it means that reaction is of the type ABSORPTIVE =
   * ALL_MOLECULES or ABSORPTIVE = ALL_VOLUME_MOLECULES since other cases
   * (REFLECTIVE/TRANSPARENT) are taken care above. But there are no products
   * for this reaction, so we do no need to go into "outcome_products()"
   * function. */

  assert(rxn_class->is_bimol());
  species_id_t all_molecules_id = p.get_all_species().get_all_molecules_species_id();
  species_id_t all_volume_molecules_id = p.get_all_species().get_all_volume_molecules_species_id();

  int result;
  bool keep_reacA = true, keep_reacB = true;

  // expecting that the surface is always the second reactant
  assert(p.get_species(rxn_class->reactant_ids[1]).is_reactive_surface());

  if (rxn_class->reactant_ids[0] == all_molecules_id ||
      rxn_class->reactant_ids[0] == all_volume_molecules_id) {
    assert(rxn_class->get_num_reactions() == 1);
    keep_reacA = false;
    result = RX_DESTROY;
  }
  else {
    // might return RX_BLOCKED
    collision.rxn_class = rxn_class;
    result = outcome_products_random(p, collision, time, pathway_index, keep_reacA, keep_reacB);
    assert(keep_reacB && "We are keeping the surface");
  }

  if (result == RX_BLOCKED) {
    return RX_A_OK; /* reflect the molecule */
  }

  if (!keep_reacA) {
#ifdef DEBUG_RXNS
    const Molecule& reacA = p.get_m(collision.diffused_molecule_id);
    DUMP_CONDITION4(
      if (!keep_reacA) {
        reacA.dump(p, "", "  defunct m:", world->get_current_iteration(), 0, false);
      }
    );
#endif

    Molecule& m = p.get_m(collision.diffused_molecule_id);
    p.set_molecule_as_defunct(m);
    return RX_DESTROY;
  }
  else {
    return result; /* RX_A_OK or RX_FLIP */
  }
}


// might return RX_BLOCKED if reaction cannot occur,
// returns 0 if positions were found
int DiffuseReactEvent::find_surf_product_positions(
    Partition& p,
    const Collision& collision,
    const BNG::RxnRule* rxn,
    const Molecule* reacA, const bool keep_reacA,
    const Molecule* reacB, const bool keep_reacB,
    const Molecule* surf_reac,
    const RxnProductsVector& actual_products,
    GridPosVector& assigned_surf_product_positions,
    uint& num_surface_products,
    bool& surf_pos_reacA_is_used
) {

  assigned_surf_product_positions.clear();

  surf_pos_reacA_is_used = keep_reacA;

  uint needed_surface_positions = 0;
  for (const ProductSpeciesIdWIndices& prod: actual_products) {
    if (p.get_species(prod.product_species_id).is_surf()) {
      needed_surface_positions++;
    }
  }
  num_surface_products = needed_surface_positions;

  assert(reacA != nullptr);
  uint num_reactants = (reacB == nullptr) ? 1 : 2;

  small_vector<GridPos> recycled_surf_prod_positions; // this array contains information on where to place the surface products
  uint initiator_recycled_index = INDEX_INVALID;

  /* list of the restricted regions for the reactants by wall */
  RegionIndicesSet rlp_wall_1, rlp_wall_2;
  /* list of the restricted regions for the reactants by object */
  RegionIndicesSet rlp_obj_1, rlp_obj_2;

  int sm_bitmask = RegionUtils::determine_molecule_region_topology(
      p, reacA, reacB, rxn->is_unimol(),
      rlp_wall_1, rlp_wall_2, rlp_obj_1, rlp_obj_2);

  // find which tiles can be recycled
  if (reacA->is_surf()) {
    if (!keep_reacA) {
      recycled_surf_prod_positions.push_back( GridPos::make_with_mol_pos(*reacA) );
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
      recycled_surf_prod_positions.push_back( GridPos::make_with_mol_pos(*reacB) );
      if (reacB->id == surf_reac->id) {
        initiator_recycled_index = recycled_surf_prod_positions.size() - 1;
      }
    }
    else if (reacB->is_surf() && keep_reacB) {
      // reacB is kept
      needed_surface_positions--;
    }
  }

  if (needed_surface_positions == 0) {
    return 0;
  }
  assigned_surf_product_positions.resize(actual_products.size());

  // do we need more tiles?
  TileNeighborVector vacant_neighbor_tiles;

  // these locations are set only for rxns with surface classes
  Vec2 hit_wall_pos2d;
  tile_index_t hit_wall_tile_index = UINT_INVALID;

  if (needed_surface_positions > recycled_surf_prod_positions.size()) {

    // find neighbors for the first surface reactant or the point where we hit the wall
    TileNeighborVector neighbor_tiles;

    if (surf_reac != nullptr) {
      // rxn with surface molecule
      Wall& wall = p.get_wall(surf_reac->s.wall_index);

      GridUtils::find_neighbor_tiles(p, surf_reac, wall, surf_reac->s.grid_tile_index, true, false, neighbor_tiles);
    }
    else {
      // rxn with surface class
      assert(collision.is_wall_collision());
      Wall& wall = p.get_wall(collision.colliding_wall_index);
      if (!wall.grid.is_initialized()) {
        wall.grid.initialize(p, wall);
      }

      hit_wall_pos2d = GeometryUtils::xyz2uv(p, collision.pos, wall);
      hit_wall_tile_index = GridUtils::uv2grid_tile_index(hit_wall_pos2d, wall);

      GridUtils::find_neighbor_tiles(p, nullptr, wall, hit_wall_tile_index, true, false, neighbor_tiles);
    }

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

  // assignment of positions
  uint num_tiles_to_recycle = min(actual_products.size(), recycled_surf_prod_positions.size());
  if (rxn->is_intermembrane_surf_rxn()) {
    // special case limited to A@C1 + B@C2 -> C@C1 + D@C2
    release_assert(needed_surface_positions == 2 && num_tiles_to_recycle == 2);
    release_assert(rxn->products.size() == 2);

    BNG::compartment_id_t compartment_prod0 = rxn->products[0].get_primary_compartment_id();

    BNG::compartment_id_t compartment_reacA = p.get_species(reacA->species_id).get_primary_compartment_id();

    // use compartment to determine location (
    if (compartment_reacA == compartment_prod0) {
      // in the same order
      assigned_surf_product_positions[0] = recycled_surf_prod_positions[0];
      assigned_surf_product_positions[0].set_reac_type(GridPosType::REACA_UV);
      assigned_surf_product_positions[1] = recycled_surf_prod_positions[1];
      assigned_surf_product_positions[1].set_reac_type(GridPosType::REACB_UV);
    }
    else {
#ifndef NDEBUG
      BNG::compartment_id_t compartment_prod1 = rxn->products[1].get_primary_compartment_id();
      assert(compartment_reacA == compartment_prod1);
#endif
      // switched order
      assigned_surf_product_positions[0] = recycled_surf_prod_positions[1];
      assigned_surf_product_positions[0].set_reac_type(GridPosType::REACB_UV);
      assigned_surf_product_positions[1] = recycled_surf_prod_positions[0];
      assigned_surf_product_positions[1].set_reac_type(GridPosType::REACA_UV);
    }
  }
  else if (needed_surface_positions == 1 && num_tiles_to_recycle == 1 && recycled_surf_prod_positions.size() >= 1) {
    if (initiator_recycled_index == INDEX_INVALID) {
      assert(recycled_surf_prod_positions.size() == 1);
      assigned_surf_product_positions[0] = recycled_surf_prod_positions[0];
    }
    else {
      assigned_surf_product_positions[0] = recycled_surf_prod_positions[initiator_recycled_index];
    }
    assigned_surf_product_positions[0].set_reac_type(GridPosType::REACA_UV);
    surf_pos_reacA_is_used = true;
  }
  else {
    uint next_available_index = 0;
    uint num_players = actual_products.size() + num_reactants;

    // assign recycled positions to products
    while (next_available_index < num_tiles_to_recycle) {
      // we must have the same number of random calls as in mcell3...
      uint rnd_num = rng_uint(&world->rng) % num_players;

      // continue until we got an index of a product
      if (rnd_num < num_reactants) {
        continue;
      }

      uint product_index = rnd_num - num_reactants;
      assert(product_index < actual_products.size());

      // we care only about surface molecules
      if (p.get_species(actual_products[product_index].product_species_id).is_vol()) {
        continue;
      }

      // skip products that we already set
      if (assigned_surf_product_positions[product_index].is_assigned()) {
        continue;
      }

      // set position for product with product_index
      assigned_surf_product_positions[product_index] = recycled_surf_prod_positions[next_available_index];
      if (assigned_surf_product_positions[product_index].has_same_wall_and_grid(*reacA)) {
        assigned_surf_product_positions[product_index].set_reac_type(GridPosType::REACA_UV);
        surf_pos_reacA_is_used = true;
      }
      else if (assigned_surf_product_positions[product_index].has_same_wall_and_grid(*reacB)) {
        assigned_surf_product_positions[product_index].set_reac_type(GridPosType::REACB_UV);
      }
      else {
        assert(false);
      }
      next_available_index++;
    }

    small_vector<bool> used_vacant_tiles;
    used_vacant_tiles.resize(vacant_neighbor_tiles.size(), false);

    // assign the first surface product of vol+wall rxn to the position where we hit the wall
    // MCell3 calls RNG and also allows to have multiple molecules to be placed on the same location,
    // MCell4 allows only a single location
    if (rxn->is_reactive_surface_rxn()) {
      // get the first unassigned surface product
      uint unassigned_product_index;
      bool found_product_index = false;
      for (unassigned_product_index = 0; unassigned_product_index < actual_products.size(); unassigned_product_index++) {
        if (rxn->products[unassigned_product_index].is_surf() &&
            !assigned_surf_product_positions[unassigned_product_index].is_assigned()) {
          found_product_index = true;
          break;
        }
      }
      release_assert(found_product_index);

      WallTileIndexPair wall_tile_indices = WallTileIndexPair(collision.colliding_wall_index, hit_wall_tile_index);
#ifndef NDEBUG
      // hit location must not be in vacant_neighbor_tiles
      uint vacant_tile_index;
      bool found_vacant_tile = false;
      for (vacant_tile_index = 0; vacant_tile_index < vacant_neighbor_tiles.size(); vacant_tile_index++) {
        if (vacant_neighbor_tiles[vacant_tile_index] == wall_tile_indices) {
          found_vacant_tile = true;
          break;
        }
      }
      assert(!found_vacant_tile);
#endif

      // call rng to have the same output as MCell3
      rng_uint(&world->rng);

      // assign position if the location is available
      const Wall& w = p.get_wall(wall_tile_indices.wall_index);
      if (w.grid.get_molecule_on_tile(hit_wall_tile_index) == MOLECULE_ID_INVALID) {
        // position was previously computed from wall and tile index
        assigned_surf_product_positions[unassigned_product_index] =
            GridPos::make_with_pos(hit_wall_pos2d, wall_tile_indices);
        assigned_surf_product_positions[unassigned_product_index].set_reac_type(GridPosType::POS_UV);
      }
    }

    // all other products are placed on one of the randomly chosen vacant tiles
    uint num_vacant_tiles = vacant_neighbor_tiles.size();
    for (uint product_index = 0; product_index < actual_products.size(); product_index++) {

      if (assigned_surf_product_positions[product_index].is_assigned()) {
        continue;
      }

      uint num_attempts = 0;
      bool found = false;
      while (!found && num_attempts < SURFACE_DIFFUSION_RETRIES) {
        assert(num_vacant_tiles != 0);
        uint rnd_num = rng_uint(&world->rng) % num_vacant_tiles;

        // is this vacant tile already used?
        if (used_vacant_tiles[rnd_num]) {
          num_attempts++;
          continue;
        }

        WallTileIndexPair grid_tile_index_pair = vacant_neighbor_tiles[rnd_num];

        // make sure we can get to the tile given the surface regions defined in the model
        if (!RegionUtils::product_tile_can_be_reached(p, grid_tile_index_pair.wall_index,
            rxn->is_unimol(), sm_bitmask, rlp_wall_1, rlp_wall_2, rlp_obj_1, rlp_obj_2)) {

          // we do not want to be checking this tile anymore
          used_vacant_tiles[rnd_num] = true;
          num_attempts++;
          continue;
        }

        assigned_surf_product_positions[product_index] = GridPos::make_without_pos(grid_tile_index_pair);
        assigned_surf_product_positions[product_index].set_reac_type(GridPosType::RANDOM);
        used_vacant_tiles[rnd_num] = true;
        found = true;
      }
      if (num_attempts >= SURFACE_DIFFUSION_RETRIES) {
        return RX_BLOCKED;
      }
    }
  }

  return 0;
}


static void update_vol_mol_after_rxn_with_surf_mol(
    Partition& p,
    const Wall& wall,
    const orientation_t product_orientation,
    const Collision& collision,
    Molecule& vm
) {
  pos_t bump = (product_orientation > 0) ? POS_EPS : -POS_EPS;
  Vec3 displacement = Vec3(2 * bump) * wall.normal;
  Vec3 new_pos_after_diffuse;

  DiffusionUtils::tiny_diffuse_3D(p, vm, displacement, wall.index, new_pos_after_diffuse);

  // update position and subpart if needed
  vm.v.pos = new_pos_after_diffuse;
  vm.v.subpart_index = p.get_subpart_index(vm.v.pos);
  p.update_molecule_reactants_map(vm);
}


void DiffuseReactEvent::handle_rxn_callback(
    Partition& p,
    const Collision& collision,
    const float_t time,
    const BNG::RxnRule* rxn,
    // reac1 is the diffused molecule and reac2 is the optional second reactant
    const Molecule* reacA,
    const Molecule* reacB,
    const MoleculeIdsVector& product_ids
) {

  // check callback
  if (world->get_callbacks().needs_rxn_callback(rxn->id)) {
    shared_ptr<API::ReactionInfo> info = make_shared<API::ReactionInfo>();

    const Molecule* reac1;
    const Molecule* reac2;

    // first reactant id that we are sending back must be always the diffusing molecule
    if (reacA->id == collision.diffused_molecule_id) {
      reac1 = reacA;
      reac2 = reacB;
    }
    else {
      assert(reacB != nullptr);
      reac1 = reacB;
      reac2 = reacA;
    }
    assert(reac2 == nullptr || reac2->id == collision.colliding_molecule_id);

    // determine type
    if (reac2 == nullptr) {
      info->type = reac1->is_vol() ?
          API::ReactionType::UNIMOL_VOLUME :
          API::ReactionType::UNIMOL_SURFACE;
    }
    else {
      if (reac1->is_vol()) {
        info->type = reac2->is_vol() ?
            API::ReactionType::VOLUME_VOLUME :
            API::ReactionType::VOLUME_SURFACE;
      }
      else {
        // surface-volume ordering or reactants is not allowed
        assert(reac2->is_surf());
        info->type = API::ReactionType::SURFACE_SURFACE;
      }
    }

    info->reactant_ids.push_back(reac1->id);
    if (reac2 != nullptr) {
      info->reactant_ids.push_back(reac2->id);
    }

    info->product_ids.insert(info->product_ids.begin(), product_ids.begin(), product_ids.end());

    info->time = time;
    info->rxn_rule_id = rxn->id;

    // pos3d
    if (reac1->is_vol()) {
      info->pos3d = collision.pos;
    }
    else {
      // collision.pos is not valid for unimol surf or surf-surf reactions
      const Wall& w = p.get_wall(reac1->s.wall_index);
      const Vec3& v0 = p.get_wall_vertex(w, 0);
      info->pos3d = GeometryUtils::uv2xyz(reac1->s.pos, w, v0);
    }

    if (rxn->is_surf_rxn()) {
      const Molecule* first_surf_reac = reac1->is_surf() ? reac1 : reac2;
      // use the first surface reactant for the first surface location
      info->pos2d = first_surf_reac->s.pos;
      info->geometry_object_id = p.get_wall(first_surf_reac->s.wall_index).object_id;
      info->partition_wall_index = first_surf_reac->s.wall_index;
    }
    else if (rxn->is_reactive_surface_rxn()) {
      const Wall& w = p.get_wall(collision.colliding_wall_index);
      info->pos2d = GeometryUtils::xyz2uv(p, collision.pos, w);
      info->geometry_object_id = w.object_id;
      info->partition_wall_index = collision.colliding_wall_index;
    }

    world->get_callbacks().do_rxn_callback(info);
  }
}


orientation_t DiffuseReactEvent::determine_orientation_depending_on_surf_comp(
    const species_id_t prod_species_id, const Molecule* surf_reac) {
  assert(surf_reac != nullptr);

  // get compartment of the surface molecule
  const Species& surf_species = world->get_all_species().get(surf_reac->species_id);
  BNG::compartment_id_t surf_comp_id = surf_species.get_primary_compartment_id();
  release_assert(surf_comp_id != BNG::COMPARTMENT_ID_NONE && "Invalid compartments used in rxn in form V(s!1).S(v!1) -> V(s) + S(v)");
  const BNG::Compartment& surf_comp = world->bng_engine.get_data().get_compartment(surf_comp_id);

  // get compartment of the product
  const Species& prod_species = world->get_all_species().get(prod_species_id);
  BNG::compartment_id_t prod_comp_id = prod_species.get_primary_compartment_id();
  release_assert(prod_comp_id != BNG::COMPARTMENT_ID_NONE && "Invalid compartments used in rxn in form V(s!1).S(v!1) -> V(s) + S(v)");

  // is the product's compartment a child or parent?
  if (surf_comp.parent_compartment_id == prod_comp_id) {
    return ORIENTATION_UP;
  }
  else if (surf_comp.children_compartments.count(prod_comp_id) != 0) {
    return ORIENTATION_DOWN;
  }
  else {
    release_assert(false && "Invalid compartments used in rxn in form V(s!1).S(v!1) -> V(s) + S(v)");
    return ORIENTATION_NOT_SET;
  }
}


// WARNING: might invalidate references
// might return RX_BLOCKED
// TODO: refactor, split into multiple functions
int DiffuseReactEvent::outcome_products_random(
    Partition& p,
    const Collision& collision,
    const float_t time,
    const rxn_class_pathway_index_t pathway_index,
    bool& keep_reacA,
    bool& keep_reacB,
    MoleculeIdsVector* optional_product_ids
) {
  assert(collision.is_mol_mol_reaction() ||
      collision.is_unimol_reaction() ||
      collision.is_wall_collision()
  );

#ifdef DEBUG_RXNS
  DUMP_CONDITION4(
    collision.dump(p, "Processing reaction:", p.stats.get_current_iteration(), time);
    cout << p.get_species(p.get_m(collision.diffused_molecule_id).species_id).name;
    if (collision.is_mol_mol_reaction()) {
      cout << " + " <<
          p.get_species(p.get_m(collision.colliding_molecule_id).species_id).name;
    }
    cout << "\nreaction_index: " << pathway_index << "\n";
    if (collision.rxn_class != nullptr) {
      collision.rxn_class->dump();
    }
    else {
      cout << "rxn_class is nullptr\n";
    }
  );
#endif

  Molecule* reacA = &p.get_m(collision.diffused_molecule_id);
  assert(reacA->is_vol() || reacA->is_surf());
  keep_reacA = false; // one product is the same as reacA
  assert(reacA != nullptr);

  Molecule* reacB = nullptr;
  keep_reacB = false; // one product is the same as reacB

  RxnClass* rxn_class = collision.rxn_class;
  assert(rxn_class != nullptr);
  const RxnRule* rxn = rxn_class->get_rxn_for_pathway(pathway_index);

  assert(rxn->reactants.size() == 1 || rxn->reactants.size() == 2);

  if (collision.is_wall_collision()) {
    keep_reacB = true;
    #ifndef NDEBUG
      // check that the second reactant is a reactive surface
      assert(rxn->reactants[1].is_simple());
      BNG::elem_mol_type_id_t mol_type_id = rxn->reactants[1].get_simple_species_mol_type_id();
      const BNG::ElemMolType& mt = p.bng_engine.get_data().get_elem_mol_type(mol_type_id);
      assert(mt.is_reactive_surface());
    #endif
  }

  assert(rxn != nullptr);

#ifdef DEBUG_CPLX_MATCHING
  cout << "Reaction to be executed:\n";
  rxn->dump(false); cout << "\n";
  rxn->dump(true); cout << "\n";
#endif

  // count this reaction if needed
  if (rxn->is_counted()) {
    assert(rxn->id !=  BNG::RXN_RULE_ID_INVALID);
    if (reacA->is_vol()) {
      p.inc_rxn_in_volume_occured_count(rxn->id, reacA->v.counted_volume_index);
    }
    else if (reacA->is_surf()) {
      p.inc_rxn_on_surface_occured_count(rxn->id, reacA->s.wall_index);
    }
  }

  Molecule* surf_reac = nullptr;
  bool reactants_swapped = false;

  // the second reactant might be a surface
  uint num_mol_reactants =
      (collision.is_wall_collision() || rxn->is_absorptive_region_rxn()) ?
      1 : rxn->reactants.size();

  if (num_mol_reactants == 2) {
    reacB = &p.get_m(collision.colliding_molecule_id);

    if (reacA->is_surf()) {
      surf_reac = reacA;
    }
    else if (reacB->is_surf()) {
      surf_reac = reacB;
    }

    // Ensure that reacA and reacB are sorted in the same order as the rxn players.
    // Needed to maintain the same behavior as in mcell3
    // With BNG, a pattern can match both of the reactants so we must check both
    // Rules are usually short so this should be rather cheap
    if (!p.bng_engine.matches_ignore_orientation(rxn->reactants[0], reacA->species_id) ||
        !p.bng_engine.matches_ignore_orientation(rxn->reactants[1], reacB->species_id)
    ) {
      Molecule* tmp_mol = reacA;
      reacA = reacB;
      reacB = tmp_mol;
      reactants_swapped = true;
    }
    // both reactants must match
    assert(p.bng_engine.matches_ignore_orientation(rxn->reactants[0], reacA->species_id));
    assert(p.bng_engine.matches_ignore_orientation(rxn->reactants[1], reacB->species_id));

    keep_reacB = rxn->is_simple_cplx_reactant_on_both_sides_of_rxn_w_identical_compartments(1);
  }
  else {
    surf_reac = reacA->is_surf() ? reacA : nullptr;
  }
  assert(p.bng_engine.matches_ignore_orientation(rxn->reactants[0], reacA->species_id));
  keep_reacA = rxn->is_simple_cplx_reactant_on_both_sides_of_rxn_w_identical_compartments(0);

  bool is_orientable = reacA->is_surf() || (reacB != nullptr && reacB->is_surf()) || collision.is_wall_collision();

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

  const RxnProductsVector& actual_products = collision.rxn_class->get_rxn_products_for_pathway(pathway_index);
  release_assert(actual_products.size() <= rxn->products.size()); // some rxn products may be still connected

  /* If the reaction involves a surface, make sure there is room for each product. */
  GridPosVector assigned_surf_product_positions; // this array contains information on where to place the surface products
  uint num_surface_products;
  bool surf_pos_reacA_is_used;
  if (is_orientable) {
    int res = find_surf_product_positions(
        p, collision, rxn, reacA, keep_reacA, reacB, keep_reacB, surf_reac, actual_products,
        assigned_surf_product_positions, num_surface_products, surf_pos_reacA_is_used);
    if (res == RX_BLOCKED) {
      return RX_BLOCKED;
    }
  }

  // free up tiles that we are probably going to reuse
  if (reacA->is_surf() && !keep_reacA) {
    p.get_wall(reacA->s.wall_index).grid.reset_molecule_tile(reacA->s.grid_tile_index);
    reacA->s.grid_tile_index = TILE_INDEX_INVALID;
  }

  if (reacB != nullptr && reacB->is_surf() && !keep_reacB) {
    p.get_wall(reacB->s.wall_index).grid.reset_molecule_tile(reacB->s.grid_tile_index);
    reacB->s.grid_tile_index = TILE_INDEX_INVALID;
  }

  // remember reactant IDs
  molecule_id_t reacA_id = reacA->id;
  molecule_id_t reacB_id = (reacB != nullptr) ? reacB->id : MOLECULE_ID_INVALID;
  molecule_id_t surf_reac_id = (surf_reac != nullptr) ? surf_reac->id : MOLECULE_ID_INVALID;

  // create and place each product
  uint current_surf_product_position_index = 0;

  int res = RX_A_OK;

  // need to determine orientations before the products are created because mcell3 does this earlier as well
  vector<orientation_t> product_orientations;
  if (is_orientable) {
    for (uint product_index = 0; product_index < rxn->products.size(); product_index++) {
      const BNG::Cplx& product = rxn->products[product_index];
      if (product.get_orientation() == ORIENTATION_NONE) {
        product_orientations.push_back( (rng_uint(&world->rng) & 1) ? ORIENTATION_UP : ORIENTATION_DOWN);
      }
      else {
        orientation_t orient = product.get_orientation();

        // set orientation relative to the durface comparmtment?
        if (orient == ORIENTATION_DEPENDS_ON_SURF_COMP) {
          orient = determine_orientation_depending_on_surf_comp(actual_products[product_index].product_species_id, surf_reac);
        }
        // flip orientation? e.g. such as for MCell3 rxn v' + s, -> s, + r, with reactants v + s'
        else if (rxn->is_bimol() && rxn->is_surf_rxn() && !rxn->is_reactive_surface_rxn()) {
          // see if orientations of the surface reactants from rule are different from the reactants
          int reacA_match = 1;
          if (reacA->is_surf() && rxn->reactants[0].get_orientation() != ORIENTATION_NONE &&
              reacA->s.orientation != rxn->reactants[0].get_orientation()) {
            reacA_match = -1;
          }
          int reacB_match = 1;
          assert(reacB != nullptr);
          if (reacB->is_surf() && rxn->reactants[1].get_orientation() != ORIENTATION_NONE &&
              reacB->s.orientation != rxn->reactants[1].get_orientation()) {
            reacB_match = -1;
          }

          // flip orientation as needed
          orient = orient * reacA_match * reacB_match;
        }

        product_orientations.push_back(orient);
      }
    }
  }
  else {
    // no surface reactant, all products will have orientataion none
    product_orientations.resize(rxn->products.size(), ORIENTATION_NONE);
  }

  MoleculeIdsVector product_ids;

  for (uint product_index = 0; product_index < actual_products.size(); product_index++) {
    const ProductSpeciesIdWIndices& actual_product = actual_products[product_index];

    // first we must check whether we are mapping a single product onto multiple complexes from the right-had side of the rule
    assert(!actual_product.rule_product_indices.empty());
    bool actual_prod_is_single_rxn_prod = actual_product.rule_product_indices.size() == 1;
    release_assert(!actual_prod_is_single_rxn_prod || *actual_product.rule_product_indices.begin() == product_index);

    // valid only if actual_prod_is_single_rxn_prod
    uint single_rxn_product_index = product_index;

    orientation_t product_orientation = ORIENTATION_NOT_SET;
    for (uint rule_product_index: actual_product.rule_product_indices) {
      const BNG::Cplx& rule_product = rxn->products[rule_product_index];
      release_assert(
          (product_orientation == ORIENTATION_NOT_SET ||
           product_orientation == product_orientations[rule_product_index]) &&
          "If a single actual product corresponds to multiple rxn products (such as in case of bond breakage), "
          "orientations of all rxn products must be the same"
      );
      product_orientation = product_orientations[rule_product_index];
    }
    release_assert(product_orientation != ORIENTATION_NOT_SET);

    // do not create anything new when the reactant is kept -
    // for bimol reactions - the diffusion simply continues
    // for unimol reactions - the unimol action action starts diffusion for the remaining timestep
    if (actual_prod_is_single_rxn_prod && rxn->is_simple_cplx_product_on_both_sides_of_rxn_w_identical_compartments(single_rxn_product_index)) {
      uint reactant_index;
      bool ok = rxn->get_assigned_cplx_reactant_for_product(single_rxn_product_index, true, reactant_index);
      assert(reactant_index == 0 || reactant_index == 1);

      Molecule* reactant = ((reactant_index == 0) ? reacA : reacB);

      // we are keeping the molecule, only the orientation changes?
      if (rxn->reactants[reactant_index].get_orientation() != product_orientation ||
          (reactant->is_surf() && reactant->s.orientation != product_orientation)) {
        // initiator volume molecule passes through wall or was flipped?
        // any surf mol -> set new orient

        // if not swapped, reacA is the diffused molecule and matches the reactant with index 0
        // reacA  matches reactant with index 0, reacB matches reactant with index 1
        assert(reactant != nullptr);
        assert(p.bng_engine.matches_ignore_orientation(rxn->products[single_rxn_product_index], reactant->species_id));

        if (reactant->is_vol()) {
          Molecule* initiator = (!reactants_swapped) ? reacA : reacB;
          if (reactant == initiator) {
            res = RX_FLIP;
          }
        }
        else if (reactant->is_surf()) {
          // surface mol - set new orientation
          reactant->s.orientation = product_orientation;
        }
        else {
          assert(false);
        }
      }

      // reactant was kept
      product_ids.push_back(reactant->id);

      continue;
    }

    species_id_t product_species_id = actual_products[product_index].product_species_id;
    const BNG::Species& species = p.get_species(product_species_id);

    molecule_id_t new_m_id;

    // set only for new vol mols when one of the reactants is surf, invalid by default
    WallTileIndexPair where_is_vm_created;

    if (species.is_vol()) {
      // create and place a volume molecule
      Vec3 pos;
      if (!collision.has_pos()) {
        // only surf-surf rxns don't have position
        assert(reacA->is_surf());
        const Wall& w_pos = p.get_wall(reacA->s.wall_index);
        pos = GeometryUtils::uv2xyz(reacA->s.pos, w_pos, p.get_wall_vertex(w_pos, 0));
      }
      else {
        pos = collision.pos;
      }
      Molecule vm_initialization(MOLECULE_ID_INVALID, product_species_id, pos, time);

      assert(is_orientable
          || (collision.type != CollisionType::VOLMOL_SURFMOL && collision.type != CollisionType::SURFMOL_SURFMOL)
      );

      Wall* wall = nullptr;
      if (is_orientable) {
        if (surf_reac != nullptr) {
          wall = &p.get_wall(surf_reac->s.wall_index);
          // position is used to schedule a diffusion action
          where_is_vm_created = surf_reac_wall_tile;
        }
        else {
          assert(collision.is_wall_collision());
          wall = &p.get_wall(collision.colliding_wall_index);

          Vec2 hit_wall_pos2d = GeometryUtils::xyz2uv(p, collision.pos, *wall);
          if (!wall->grid.is_initialized()) {
            wall->initialize_grid(p);
          }
          tile_index_t hit_wall_tile_index = GridUtils::uv2grid_tile_index(hit_wall_pos2d, *wall);

          where_is_vm_created = WallTileIndexPair(collision.colliding_wall_index, hit_wall_tile_index);
        }

        // tiny diffuse done in update_vol_mol_after_rxn_with_surf_mol
        // cannot cross walls
        CollisionUtils::update_counted_volume_id_when_crossing_wall(
            p, *wall, product_orientation, vm_initialization);
      }

      // adding molecule might invalidate references of already existing molecules and also of species
      Molecule& new_vm = p.add_volume_molecule(vm_initialization);

      // id used to schedule a diffusion action
      new_m_id = new_vm.id;

      const BNG::Species& species_new_ref = p.get_species(product_species_id);
      new_vm.set_flag(MOLECULE_FLAG_VOL);
      new_vm.set_flag(MOLECULE_FLAG_SCHEDULE_UNIMOL_RXN);

      if (is_orientable) {
        // - for an orientable reaction, we need to move products away from the surface
        //   to ensure they end up on the correct side of the plane.
        update_vol_mol_after_rxn_with_surf_mol(
            p, *wall, product_orientation, collision, new_vm
        );
      }
    #ifdef DEBUG_RXNS
      DUMP_CONDITION4(
        new_vm.dump(p, "", "  created vm:", world->get_current_iteration(), time);
      );
    #endif
    }
    else {
      // see release_event_t::place_single_molecule_onto_grid, merge somehow

      // get info on where to place the product and increment the counter
      assert(current_surf_product_position_index < assigned_surf_product_positions.size());
      const GridPos& new_grid_pos = assigned_surf_product_positions[current_surf_product_position_index];
      assert(new_grid_pos.is_assigned());

      Vec2 pos;
      switch (new_grid_pos.type) {
        case GridPosType::REACA_UV:
          if (rxn->is_unimol() && (num_surface_products == 2) && (surf_reac != NULL)) {
            // there must be a single position with random setting for the second surface product,
            // find and use it
            const GridPos* second_prod_grid_pos =
                GridPos::get_second_surf_product_pos(assigned_surf_product_positions, current_surf_product_position_index);
            assert(second_prod_grid_pos != nullptr);

            pos = GridPosition::find_closest_position(p, new_grid_pos, *second_prod_grid_pos);
          }
          else {
            assert(new_grid_pos.pos_is_set);
            pos = new_grid_pos.pos;
          }
          break;
        case GridPosType::REACB_UV:
        case GridPosType::POS_UV:
          assert(new_grid_pos.pos_is_set);
          pos = new_grid_pos.pos;
          break;
        case GridPosType::RANDOM:
          if (rxn->is_unimol() && !surf_pos_reacA_is_used && (num_surface_products == 2)) {
            const GridPos* second_prod_grid_pos =
                GridPos::get_second_surf_product_pos(assigned_surf_product_positions, current_surf_product_position_index);
            assert(second_prod_grid_pos != nullptr);

            pos = GridPosition::find_closest_position(p, new_grid_pos, *second_prod_grid_pos);
          }
          else {
            const Wall& wall = p.get_wall(new_grid_pos.wall_index);
            pos = GridUtils::grid2uv_random(wall, new_grid_pos.tile_index, world->rng);
          }
          break;
        default:
          release_assert(false);
      }

      // create our new molecule
      Molecule sm_to_add(MOLECULE_ID_INVALID, product_species_id, pos, time);
      sm_to_add.set_flag(MOLECULE_FLAG_SURF);
      sm_to_add.set_flag(MOLECULE_FLAG_SCHEDULE_UNIMOL_RXN);

      // set wall and grid information
      sm_to_add.s.wall_index = new_grid_pos.wall_index;
      sm_to_add.s.grid_tile_index = new_grid_pos.tile_index;

      // and orientation
      sm_to_add.s.orientation = product_orientation;

      // might invalidate references of already existing molecules
      // returned molecule has its id set
      Molecule& new_sm = p.add_surface_molecule(sm_to_add);
      Grid& grid = p.get_wall(new_grid_pos.wall_index).grid;
      grid.set_molecule_tile(new_sm.s.grid_tile_index, new_sm.id);
      new_m_id = new_sm.id;

      #ifdef DEBUG_RXNS
        DUMP_CONDITION4(
          new_sm.dump(p, "", "  created sm:", world->get_current_iteration(), time);
        );
      #endif

      current_surf_product_position_index++;
    }

    p.get_m(new_m_id).diffusion_time = time;
    if (before_this_iterations_end(time)) {
      new_diffuse_actions.push_back(DiffuseAction(new_m_id, where_is_vm_created));
    }

    product_ids.push_back(new_m_id);

    // refresh reacA and reacB pointers, we are added new molecules in this loop and they point to that vector
    reacA = &p.get_m(reacA_id);
    reacB = (reacB_id != MOLECULE_ID_INVALID) ? &p.get_m(reacB_id) : nullptr;
    surf_reac = (surf_reac_id != MOLECULE_ID_INVALID) ? &p.get_m(surf_reac_id) : nullptr;
  } // end for - product creation

  // this is a safe point where we can manipulate contents of this event
  handle_rxn_callback(p, collision, time, rxn, reacA, reacB, product_ids);

  if (optional_product_ids != nullptr) {
    *optional_product_ids = product_ids;
  }

  // we might need to swap info on which reactant was kept
  if (reactants_swapped) {
    bool tmp = keep_reacA;
    keep_reacA = keep_reacB;
    keep_reacB = tmp;
  }

  return res;
}

// ---------------------------------- unimolecular reactions ----------------------------------

// WARNING: might invalidate molecule and species references
// returns true if molecule survived
bool DiffuseReactEvent::outcome_unimolecular(
    Partition& p,
    Molecule& m,
    const float_t scheduled_time,
    RxnClass* rxn_class,
    const rxn_class_pathway_index_t pathway_index,
    MoleculeIdsVector* optional_product_ids
) {
  molecule_id_t id = m.id;

  // a PATHWAY_INDEX_NO_RXN pathway might have been selected
  // when no unimol rxn matched compartments
  if (pathway_index >= PATHWAY_INDEX_LEAST_VALID) {

    Vec3 pos;
    if (m.is_vol()) {
      pos = m.v.pos;
    }
    else if (m.is_surf()) {
      Wall& w = p.get_wall(m.s.wall_index);
      const Vec3& wall_vert0 = p.get_geometry_vertex(w.vertex_indices[0]);
      pos = GeometryUtils::uv2xyz(m.s.pos, w, wall_vert0);
    }
    else {
      pos = Vec3(POS_INVALID);
      assert(false);
    }

    Collision collision(
        CollisionType::UNIMOLECULAR, &p, m.id, scheduled_time, pos, rxn_class);

    bool ignoredA, ignoredB;
    // creates new molecule(s) as output of the unimolecular reaction
    // !! might invalidate references (we might reorder defuncting and outcome call later)
    int outcome_res = outcome_products_random(
        p, collision, scheduled_time, pathway_index, ignoredA, ignoredB, optional_product_ids);
    assert(outcome_res == RX_A_OK || outcome_res == RX_BLOCKED);

    Molecule& m_new_ref = p.get_m(id);

    const RxnRule* unimol_rx = rxn_class->get_rxn_for_pathway(pathway_index);

    // and defunct this molecule if it was not kept
    assert(unimol_rx->reactants.size() == 1 || unimol_rx->is_absorptive_region_rxn());
    if (outcome_res != RX_BLOCKED && !unimol_rx->is_simple_cplx_reactant_on_both_sides_of_rxn_w_identical_compartments(0)) {
#ifdef DEBUG_RXNS
      DUMP_CONDITION4(
        m_new_ref.dump(p, "", m_new_ref.is_vol() ? "Unimolecular vm defunct:" : "Unimolecular sm defunct:", world->get_current_iteration(), scheduled_time, false);
      );
#endif
      p.set_molecule_as_defunct(m_new_ref);
      return false;
    }
  }

  // molecule survived
  // we must reschedule the molecule's unimol rxn, this will happen right away
  // during the molecule's diffusion
  Molecule& m_new_ref = p.get_m(id);
  m_new_ref.set_flag(MOLECULE_FLAG_SCHEDULE_UNIMOL_RXN);

  return true;
}


// returns true if molecule survived
bool DiffuseReactEvent::cross_transparent_wall(
    Partition& p,
    const Collision& collision,
    Molecule& vm, // moves vm to the reflection point
    Vec3& remaining_displacement,
    float_t& t_steps,
    float_t& elapsed_molecule_time,
    wall_index_t& last_hit_wall_index
) {
#ifdef DEBUG_COUNTED_VOLUMES
  vm.dump(p, "- Before crossing: ");
#endif

  const Wall& w = p.get_wall(collision.colliding_wall_index);

#ifdef DEBUG_TRANSPARENT_SURFACES
  std::cout << "Crossed a transparent wall, side: " << w.side << "\n";
#endif

  // check if we are not crossing a compartment boundary
  const GeometryObject& obj = p.get_geometry_object(w.object_index);
  if (obj.represents_compartment()) {
    molecule_id_t orig_vm_id = vm.id;

    // get species, pretty inefficient, may need to be cached if this is a feature that is used often
    BNG::Species new_species = p.get_species(vm.species_id);
    new_species.set_compartment_id(BNG::COMPARTMENT_ID_NONE);
    new_species.finalize_species(world->bng_engine.get_config(), true);
    species_id_t new_species_id = p.get_all_species().find_or_add(new_species, true);

    // using the same computation of collision time as in collide_and_react_with_surf_mol
    float_t collision_time = elapsed_molecule_time + t_steps * collision.time;

    // we create a new molecule on the boundary, move it a tiny bit in the right direction
    Molecule vm_initialization(
        MOLECULE_ID_INVALID, new_species_id,
        collision.pos + remaining_displacement * Vec3(EPS),
        event_time + collision_time);
    vm_initialization.v.previous_wall_index = w.index;
    Molecule& new_vm = p.add_volume_molecule(vm_initialization);
    new_vm.set_flag(MOLECULE_FLAG_VOL);
    new_vm.set_flag(MOLECULE_FLAG_SCHEDULE_UNIMOL_RXN);

    // schedule a diffusion action
    new_vm.diffusion_time = new_vm.birthday;
    if (before_this_iterations_end(new_vm.birthday)) {
      new_diffuse_actions.push_back(DiffuseAction(new_vm.id));
    }

    // and destroy the previous one
    Molecule& orig_vm_new_ref = p.get_m(orig_vm_id);
    p.set_molecule_as_defunct(orig_vm_new_ref);

#ifdef DEBUG_COUNTED_VOLUMES
    vm_initialization.dump(p, "- After crossing: ");
#endif
    return false;
  }
  else {
    // keep the current molecule and just move it

    // Update molecule location to the point of wall crossing
    vm.v.pos = collision.pos;
    vm.v.subpart_index = p.get_subpart_index(vm.v.pos);

    const BNG::Species& sp = p.bng_engine.get_all_species().get(vm.species_id);
    CollisionUtils::update_counted_volume_id_when_crossing_wall(
        p, w, collision.get_orientation_against_wall(), vm);

    // ignore the collision time, it is a bit earlier and does not fit for multiple collisions in the
    // same time step
    float_t t_smash = collision.time;

    remaining_displacement = remaining_displacement * Vec3(1.0 - t_smash);
    elapsed_molecule_time += t_steps * t_smash;

    t_steps *= (1.0 - t_smash);
    if (t_steps < EPS) {
      t_steps = EPS;
    }

    last_hit_wall_index = w.index;
#ifdef DEBUG_COUNTED_VOLUMES
  vm.dump(p, "- After crossing: ");
#endif

    return true;
  }
}

// ---------------------------------- dumping methods ----------------------------------

void DiffuseReactEvent::dump(const std::string ind) const {
  cout << ind << "Diffuse-react event:\n";
  std::string ind2 = ind + "  ";
  BaseEvent::dump(ind2);
  cout << ind2 << "barrier_time_from_event_time: \t\t" << time_up_to_next_barrier << " [float_t] (may be unset initally)\t\t\n";
}


} /* namespace mcell */
