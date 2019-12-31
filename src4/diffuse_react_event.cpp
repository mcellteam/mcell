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

// include implementations of utility functions
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
    diffuse_single_molecule(p, id, diffusion_time_step, WallTileIndexPair());
  }


  // 2) simulate remaining time of molecules created with reactions or
  // scheduled unimolecular reactions
  // need to call .size() each iteration because the size can increase,
  // again, we are using it as a queue and we do not follow the time when
  // they were created
  for (uint i = 0; i < new_diffuse_or_unimol_react_actions.size(); i++) {
    const DiffuseOrUnimolReactionAction& action = new_diffuse_or_unimol_react_actions[i];

    if (action.type == DiffuseOrUnimolReactionAction::Type::DIFFUSE) {
      diffuse_single_molecule(
          p, action.id,
          event_time + diffusion_time_step - action.scheduled_time,
          action.where_created_this_iteration // making a copy of this pair
      );
    }
    else {
      react_unimol_single_molecule(p, action.id, action.scheduled_time, action.unimol_rx);
    }
  }

  new_diffuse_or_unimol_react_actions.clear();
}


void DiffuseReactEvent::diffuse_single_molecule(
    Partition& p,
    const molecule_id_t m_id,
    const float_t time_up_to_event_end, // how much diffusion time we should simulate
    WallTileIndexPair wall_tile_pair_where_created_this_iteration // set only for newly created molecules
) {
  float_t event_time_end = event_time + diffusion_time_step;
  Molecule& m = p.get_m(m_id);

  if (m.is_defunct())
    return;

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
    DiffuseOrUnimolReactionAction unimol_react_action(
        DiffuseOrUnimolReactionAction::Type::UNIMOL_REACT, m.id, m.unimol_rx_time, m.unimol_rx);
    // handle this iteration
    new_diffuse_or_unimol_react_actions.push_back(unimol_react_action);
  }
  /*else {
    p.add_unimolecular_action(diffusion_time_step, unimol_react_action);
  }*/


#ifdef DEBUG_DIFFUSION
  const Species& debug_species = world->get_species(m.species_id);
  DUMP_CONDITION4(
    // the subtraction of diffusion_time_step doesn't make much sense but is needed to make the dump the same as in mcell3
    // need to check it further
    const char* title =
        (debug_species.can_diffuse()) ?
            (m.is_vol() ? "Diffusing vm:" : "Diffusing sm:") :
            (m.is_vol() ? "Not diffusing vm:" : "Not diffusing sm:");
    m.dump(p, "", title, world->get_current_iteration(), event_time_end - time_up_to_event_end);
  );
#endif

  // we might need to adjust remaining time step if this molecule has a unimolecular reaction
  // within this event's time step range
  float_t remaining_time_step;
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

  bool was_defunct;
  vec3_t new_vpos;
  vec2_t new_spos;
  subpart_index_t new_subpart_index;
  wall_index_t new_wall_index;
  if (m.is_vol()) {
    diffuse_vol_molecule(
        p, m_id, remaining_time_step,
        was_defunct, new_vpos, new_subpart_index,
        wall_tile_pair_where_created_this_iteration
    );
  }
  else {
    float_t advance_time;
    diffuse_surf_molecule(
        p, m_id, event_time_end - time_up_to_event_end, remaining_time_step,
        was_defunct, new_spos, new_wall_index,
        advance_time
    );
  }

  if (!was_defunct) { // FIXME: unify - we are changing the partition, but not the wall here
    // need to get a new reference
    Molecule& m_new_ref = p.get_m(m_id);

    if (m.is_vol()) {
      // finally move molecule to its destination
      m_new_ref.v.pos = new_vpos;

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
    const float_t remaining_time_step,
    bool& was_defunct,
    vec3_t& new_pos,
    subpart_index_t& new_subpart_index,
    WallTileIndexPair& wall_tile_pair_where_created_this_iteration
) {
  Molecule& m = p.get_m(vm_id);
  const Species& species = world->get_species(m.species_id);

  // diffuse each molecule - get information on position change
  vec3_t displacement;
  float_t r_rate_factor;
  diffusion_util::compute_vol_displacement(
      species, remaining_time_step, world->rng, displacement, r_rate_factor);

#ifdef DEBUG_DIFFUSION
  DUMP_CONDITION4(
      if (species.can_diffuse()) {
        displacement.dump("  displacement:", "");
      }
  );
#endif
  // note: we are ignoring use_expanded_list setting compared to mcell3

  // detect collisions with other molecules
  vec3_t remaining_displacement = displacement;

  RayTraceState state;
  collision_vector_t molecule_collisions;
  was_defunct = false;
  wall_index_t reflected_wall_index = WALL_INDEX_INVALID;
  float_t updated_remaining_time_step = remaining_time_step; // == t_steps
  float_t elapsed_molecule_time = diffusion_time_step - remaining_time_step; // == vm->t
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

      if (collision.is_vol_mol_collision()) {
        // ignoring immediate collisions
        if (collision.time < EPS) {
          continue;
        }

        // evaluate reaction associated with this collision
        // for now, do the change right away, but we might need to cache these changes and
        // do them after all diffusions were finished
        // warning: might invalidate references to p.volume_molecules array! returns true in that case
        if (collide_and_react_with_vol_mol(
              p, collision, remaining_displacement,
              updated_remaining_time_step, r_rate_factor)
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
          info.pos_before_hit = vm_new_ref.v.pos;

          wall_hit_callback(info, world->wall_hit_callback_clientdata);
        }
#ifdef DEBUG_WALL_COLLISIONS
        cout << "Wall collision: \n";
        const GeometryObject* geom_obj = p.get_geometry_object_if_exists(colliding_wall.object_id);
        assert(geom_obj != nullptr);
        cout << "  mol id: " << vm_new_ref.id << ", species: " << world->get_species(vm_new_ref.species_id).name << "\n";
        cout << "  obj: " << geom_obj->name << ", id: " << geom_obj->id << "\n";
        cout << "  wall id: " << colliding_wall.id << "\n";
        cout << "  time: " << float_t(world->get_current_iteration()) + collision.time << "\n";
        cout << "  pos: " << collision.pos << "\n";
#endif

        // check possible reaction with surface molecules
        if (species.has_flag(SPECIES_FLAG_CAN_VOLSURF) && colliding_wall.has_initialized_grid()) {
          int collide_res = collide_and_react_with_surf_mol(
              p, collision, updated_remaining_time_step,
              r_rate_factor, elapsed_molecule_time,
              wall_tile_pair_where_created_this_iteration
          );

          if (collide_res == 1) { // FIXME: use enum
            was_defunct = true;
            return;
          }
        }

        if (!was_defunct) {
          elapsed_molecule_time += updated_remaining_time_step * collision.time;
          // if a molecule was reflected, changes its position to the reflection point
          int res = CollisionUtil::reflect_or_periodic_bc(
              p, collision,
              vm_new_ref, remaining_displacement, updated_remaining_time_step, reflected_wall_index
          );
          assert(res == 0 && "Periodic box BCs are not supported yet");
        }

        // molecule could have been moved
        subpart_index_t subpart_after_wall_hit = p.get_subpartition_index(vm_new_ref.v.pos);
        // change subpartition if needed
        p.change_molecule_subpartition(vm_new_ref, subpart_after_wall_hit);

        break; // we reflected and did not react, do ray_trace again
      }
    }

  } while (unlikely(state != RayTraceState::FINISHED && !was_defunct));
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
  p.get_simulation_stats().inc_ray_voxel_tests();

  RayTraceState res_state = RayTraceState::FINISHED;
  collisions.clear();

  const WorldConstants& world_constants = p.get_world_constants();
  float_t radius = world_constants.rx_radius_3d;

  // first get what subpartitions might be relevant
  SubpartIndicesVector crossed_subparts_for_walls;
  subpart_indices_set_t crossed_subparts_for_molecules;
  subpart_index_t last_subpartition_index;
  CollisionUtil::collect_crossed_subparts(
      p, vm, remaining_displacement,
      radius, world_constants.subpartition_edge_length,
      true, crossed_subparts_for_walls,
      crossed_subparts_for_molecules, last_subpartition_index
  );


  vec3_t& displacement_up_to_wall_collision = remaining_displacement;

  // check wall collisions in the crossed subparitions,
  // stop at first crossing because crossed_subparts_for_walls are ordered
  // and we are sure that if we hit a wall in the actual supartition, we cannot
  // possibly hit another wall in a subparition that follows
  if (!crossed_subparts_for_walls.empty()) {
    for (subpart_index_t subpart_index: crossed_subparts_for_walls) {

      CollisionUtil::collect_wall_collisions( // mcell3 does this only for the current subvol
          p,
          vm,
          subpart_index,
          previous_reflected_wall,
          rng,
          displacement_up_to_wall_collision,
          collisions
      );

      if (!collisions.empty()) {
        res_state = RayTraceState::RAY_TRACE_HIT_WALL;
        break;
      }
    }
  }

  if (res_state == RayTraceState::RAY_TRACE_HIT_WALL) {
    // recompute collect_crossed_subparts if there was a wall collision
    crossed_subparts_for_molecules.clear();
    CollisionUtil::collect_crossed_subparts(
        p, vm, displacement_up_to_wall_collision,
        radius,
        world_constants.subpartition_edge_length,
        false, crossed_subparts_for_walls, // not filled this time
        crossed_subparts_for_molecules, last_subpartition_index
    );
  }

  // check molecule collisions for each SP
  for (subpart_index_t subpart_index: crossed_subparts_for_molecules) {
    // get cached reacting molecules for this SP
    UintSet& sp_reactants = p.get_volume_molecule_reactants(subpart_index, vm.species_id);

    // for each molecule in this SP
    for (molecule_id_t colliding_vm_id: sp_reactants) {
      CollisionUtil::collide_mol_loop_body(
          world_constants,
          p,
          vm,
          colliding_vm_id,
          displacement_up_to_wall_collision,
          radius,
          collisions
      );
    }
  }

  // the value is valid only when RAY_TRACE_FINISHED is returned
  new_subpart_index = last_subpartition_index; // FIXME: this is valid only when there was no collision
  new_pos = vm.v.pos + remaining_displacement; // FIXME: this is overwritten in case of a collision, should I do something about it?

  return res_state; // no wall was hit
}


// handle collision of two volume molecules: checks probability of reaction,
// executes this reaction, removes reactants and creates products
// returns true if reaction has occured and the first reactant was destroyed
bool DiffuseReactEvent::collide_and_react_with_vol_mol(
    Partition& p,
    const Collision& collision,
    vec3_t& displacement,
    const float_t remaining_time_step,
    const float_t r_rate_factor
)  {

  Molecule& colliding_molecule = p.get_m(collision.colliding_molecule_id); // am
  Molecule& diffused_molecule = p.get_m(collision.diffused_molecule_id); // m

  // returns 1 when there are no walls at all
  float_t factor = ExactDiskUtil::exact_disk(
      p, collision.pos, displacement, p.get_world_constants().rx_radius_3d,
      diffused_molecule, colliding_molecule,
      p.get_world_constants().use_expanded_list
  );

  if (factor < 0) { /* Probably hit a wall, might have run out of memory */
    return 0; /* Reaction blocked by a wall */
  }

  const Reaction& rx = *collision.rx;
  //  rx->prob_t is always NULL in out case update_probs(world, rx, m->t);
  // returns which reaction pathway to take
  float_t scaling = factor * r_rate_factor;
  int i = rx_util::test_bimolecular(
    rx, world->rng, colliding_molecule, diffused_molecule, scaling, 0);

  if (i < RX_LEAST_VALID_PATHWAY) {
    return false;
  }
  else {
    // might invalidate references
    int j = outcome_bimolecular(p, collision, i, remaining_time_step);
    assert(j == RX_DESTROY);
    return true;
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

  ReactionsVector matching_rxns;
  rx_util::trigger_bimolecular(
    *p.get_world_constants().bimolecular_reactions_map,
    diffused_molecule, colliding_molecule,
    collision_orientation, colliding_molecule.s.orientation,
    matching_rxns
  );

  if (matching_rxns.empty()) {
    return -1;
  }

  // FIXME: this code is very similar to code in react_2D_neighbors
  small_vector<float_t> scaling_coefs;
  for (size_t i = 0; i < matching_rxns.size(); i++) {
    const Reaction* rxn = matching_rxns[i];
    assert(rxn != nullptr);

    scaling_coefs.push_back(r_rate_factor / grid.binding_factor);
  }

  int selected_rx_pathway;
  int reactant_index;
  if (matching_rxns.size() == 1) {
    selected_rx_pathway = rx_util::test_bimolecular(
        *matching_rxns[0], world->rng,
        diffused_molecule, colliding_molecule,
        scaling_coefs[0], 0);

    assert(selected_rx_pathway <= 0 && "Only one pathway supported for now (with index 0)");
    reactant_index = 0;
  }
  else {
    bool all_neighbors_flag = true;
    reactant_index = rx_util::test_many_bimolecular(matching_rxns, scaling_coefs, 0, world->rng, false);
    selected_rx_pathway = 0; // TODO_PATHWAYS: use value from test_many_bimolecular
  }


  if (reactant_index == RX_NO_RX || selected_rx_pathway < RX_LEAST_VALID_PATHWAY) {
    return -1; /* No reaction */
  }

  /* run the reaction */
  Collision rx_collision = Collision(
      CollisionType::VOLMOL_SURFMOL,
      &p,
      collision.diffused_molecule_id,
      elapsed_molecule_time + remaining_time_step * collision.time,
      collision.pos,
      colliding_molecule.id,
      matching_rxns[selected_rx_pathway] // FIXME_PATHWAYS: this should contain all the pathways
  );

  int outcome_bimol_result = outcome_bimolecular(
      p, rx_collision, selected_rx_pathway, remaining_time_step
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
    const float_t current_time,
    const float_t remaining_time_step,
    bool& was_defunct,
    vec2_t& new_loc,
    wall_index_t& new_wall_index,
    float_t& advance_time
) {
  Molecule& sm = p.get_m(sm_id);
  const Species& species = world->get_species(sm.species_id);

  float_t steps = 0.0;
  float_t t_steps = 0.0;
  float_t space_factor = 0.0;


  wall_index_t original_wall_index = sm.s.wall_index;

  /* Where are we going? */
  if (species.time_step > remaining_time_step) {
    t_steps = remaining_time_step;
    steps = remaining_time_step / species.time_step ;
  }
  else {
    t_steps = species.time_step;
    steps = 1.0;
  }
  if (steps < EPS) {
    t_steps = EPS * species.time_step;
    steps = EPS;
  }

  if (species.space_step != 0) {

    if (steps == 1.0) {
      space_factor = species.space_step;
    }
    else {
      space_factor = species.space_step * sqrt_f(steps);
    }

    for (int find_new_position = (SURFACE_DIFFUSION_RETRIES + 1);
         find_new_position > 0; find_new_position--) {

      vec2_t displacement;
      diffusion_util::compute_surf_displacement(species, space_factor, world->rng, displacement);


  #ifdef DEBUG_DIFFUSION
    DUMP_CONDITION4(
        if (species.can_diffuse()) {
          displacement.dump("  displacement:", "");
        }
    );
  #endif

      assert(!species.has_flag(SPECIES_FLAG_SET_MAX_STEP_LENGTH) && "not supported yet");

      // ray_trace does the movement and all other stuff
      new_wall_index =
          ray_trace_surf(p, species, sm, displacement, new_loc, was_defunct/*, &rxp, &hd_info*/);

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
        if (diffusion_util::move_sm_on_same_triangle(p, sm, new_loc)) {
          continue;
        }
      }
      // After diffusing, we ended up on a NEW triangle.
      else {
        if (diffusion_util::move_sm_to_new_triangle(p, sm, new_loc, new_wall_index)) {
          continue;
        }
      }
      find_new_position = 0;
    }
  } // if (species.space_step != 0)

  advance_time = t_steps;


  // NOTE: what about molecules that cannot diffuse?
  assert(!species.has_flag(SPECIES_FLAG_CAN_SURFSURFSURF) && "Not supported");
  if (species.has_flag(SPECIES_FLAG_CAN_SURFSURF)) {
    assert(!species.has_flag(SPECIES_FLAG_CANT_INITIATE) && "Not sure what to do here");

    if (react_2D_all_neighbors(p, sm, current_time, advance_time)) {
      was_defunct = true;
    }
  }


  // reactions in react_2D_all_neighbors could have invalidated
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


// returns true if molecule still exists
bool DiffuseReactEvent::react_2D_all_neighbors(
    Partition& p,
    Molecule& sm,
    const float_t current_time,
    const float_t remaining_time_step
) {

  const Wall& wall = p.get_wall(sm.s.wall_index);

  TileNeighborVector neighbors;
  GridUtil::find_neighbor_tiles(p, sm, wall, sm.s.grid_tile_index, false, /* true,*/ neighbors);

  if (neighbors.empty()) {
    return true;
  }

  const Species& sm_species = world->get_species(sm.species_id);


  size_t l = 0;
  // array, each item corresponds to one potential reaction
  small_vector<float_t> correction_factors;
  small_vector<molecule_id_t> reactant_molecule_ids;
  ReactionsVector matching_rxns;

  /* step through the neighbors */
  for (const WallTileIndexPair& neighbor: neighbors) {
    Grid& ngrid = p.get_wall(neighbor.wall_index).grid;

    // is there something on the tile?
    // TODO_LATER: this filtering should be done already while looking for neighbors
    molecule_id_t nid = ngrid.get_molecule_on_tile(neighbor.tile_index);
    if (nid == MOLECULE_ID_INVALID) {
      continue;
    }

    Molecule& nsm = p.get_m(nid);
    const Species& nsm_species = world->get_species(nsm.species_id);

#ifdef DEBUG_REACTIONS
    DUMP_CONDITION4(
      // the subtraction of diffusion_time_step doesn't make much sense but is needed to make the dump the same as in mcell3
      // need to check it further
      nsm.dump(p, "", "  checking in react_2D_all_neighbors: ", world->get_current_iteration(), 0.0/*event_time_end - time_up_to_event_end - diffusion_time_step*//*time ???*/);
    );
#endif

    /* check whether the neighbor molecule is behind
       the restrictive region boundary   */
    assert(!sm_species.has_flag(SPECIES_FLAG_CAN_REGION_BORDER) && "TODO_LATER");
    assert(!nsm_species.has_flag(SPECIES_FLAG_CAN_REGION_BORDER) && "TODO_LATER");
    assert(!nsm_species.has_flag(SPECIES_FLAG_EXTERNAL_SPECIES) && "TODO_LATER");

    // returns value >=1 if there can be a reaction
    size_t orig_num_rxsn = matching_rxns.size();
    rx_util::trigger_bimolecular_orientation_from_mols(
        world->bimolecular_reactions_map,
        sm, nsm,
        matching_rxns
    );

    // extend arrays holding additional information
    // FIXME: the same code is in collide_and_react_with_surf_mol
    for (size_t i = orig_num_rxsn; i < matching_rxns.size(); i++) {
      const Reaction* rxn = matching_rxns[i];
      assert(rxn != nullptr);

      correction_factors.push_back(remaining_time_step / ngrid.binding_factor);
      reactant_molecule_ids.push_back(nsm.id);
    }
  }

  size_t num_matching_rxns = matching_rxns.size();
  if (num_matching_rxns == 0) {
    return true;
  }

  int selected_rx_pathway;
  Collision collision;
  float_t local_prob_factor = 3.0 / neighbors.size();
  int reactant_index;
  if (num_matching_rxns == 1) {
    // figure out what should happen
    /* Calculate local_prob_factor for the reaction probability.
       Here we convert from 3 neighbor tiles (upper probability
       limit) to the real "num_nbrs" neighbor tiles. */

    selected_rx_pathway = rx_util::test_bimolecular(
        *matching_rxns[0], world->rng,
        sm, p.get_m(reactant_molecule_ids[0]),
        correction_factors[0], local_prob_factor);

    assert(selected_rx_pathway <= 0 && "Only one pathway supported for now (with index 0)");
    reactant_index = 0;
  }
  else {
    bool all_neighbors_flag = true;
    reactant_index = rx_util::test_many_bimolecular(matching_rxns, correction_factors, local_prob_factor, world->rng, all_neighbors_flag);
    selected_rx_pathway = 0; // TODO_PATHWAYS: use value from test_many_bimolecular
  }

  if (reactant_index == RX_NO_RX || selected_rx_pathway < RX_LEAST_VALID_PATHWAY) {
    return true; /* No reaction */
  }

  collision = Collision(CollisionType::SURFMOL_SURFMOL, &p, sm.id, current_time - event_time, reactant_molecule_ids[reactant_index], matching_rxns[reactant_index]);

  /* run the reaction */
  int outcome_bimol_result = outcome_bimolecular(
      p, collision, selected_rx_pathway, remaining_time_step
  );


  return true;
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
    vec2_t& new_pos,
    bool& was_defunct
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
    // FIXME: enum for index_edge_was_hit?
    int index_edge_was_hit =
        GeometryUtil::find_edge_point(*this_wall, this_pos, this_disp, boundary_pos);

    // Ambiguous edge collision. Give up and try again from diffuse_2D.
    if (index_edge_was_hit == -2) {
      sm.s.pos = orig_pos;
      // hit_data_info = hit_data_head;
      return WALL_INDEX_INVALID;
    }

    // We didn't hit the edge. Stay inside this wall. We're done!
    else if (index_edge_was_hit == -1) {
      new_pos = this_pos + this_disp;

      // ???
      sm.s.pos = orig_pos;
      // *hit_data_info = hit_data_head;
      return this_wall->index;
    }


    // Neither ambiguous (-2) nor inside wall (-1), must have hit edge (0, 1, 2)
    vec2_t old_pos = this_pos;

    /* We hit the edge - check for the reflection/absorption from the
       edges of the wall if they are region borders
       Note - here we test for potential collisions with the region
       border while moving INSIDE OUT */
    assert(!species.has_flag(SPECIES_FLAG_CAN_REGION_BORDER) && "not supported yet");

    /* no reflection - keep going */
    vec2_t new_disp;
    wall_index_t target_wall_index =
        GeometryUtil::traverse_surface(*this_wall, old_pos, index_edge_was_hit, this_pos);

    if (target_wall_index != WALL_INDEX_INVALID) {
      assert(!species.has_flag(SPECIES_FLAG_CAN_REGION_BORDER) && "not supported yet");

      this_disp = old_pos + this_disp;

      #ifndef NDEBUG
        Edge& e = const_cast<Edge&>(this_wall->edges[index_edge_was_hit]);
        assert(e.is_initialized());
        e.debug_check_values_are_uptodate(p);
      #endif

      GeometryUtil::traverse_surface(*this_wall, this_disp, index_edge_was_hit, new_disp);
      this_disp = new_disp - this_pos;
      this_wall = &p.get_wall(target_wall_index);
      continue;
    }

    /* If we reach this point, assume we reflect off the edge since there is no
     * neighboring wall
     *
     * NOTE: this_pos has been corrupted by traverse_surface; use old_pos to find
     * out whether the present wall edge is a region border
     */
    new_disp = this_disp - (boundary_pos - old_pos);

    switch (index_edge_was_hit) {
      case 0:
        new_disp.v *= -1.0;
        break;
      case 1: {
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
      case 2: {
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
        UNHANDLED_CASE(index_edge_was_hit);
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

  const Reaction* rx = rx_util::pick_unimol_rx(world, m.species_id);
  if (rx == nullptr) {
    return;
  }

  float_t time_from_now = rx_util::compute_unimol_lifetime(p, world->rng, m, rx);

  float_t scheduled_time = curr_time + time_from_now;

  // we need to store the end time to the molecule because oit is needed in diffusion to
  // figure out whether we should do the whole time step
  m.unimol_rx_time = scheduled_time;
  m.unimol_rx = rx;
}


// based on mcell3's check_for_unimolecular_reaction
// might invalidate vm references
void DiffuseReactEvent::react_unimol_single_molecule(
    Partition& p,
    const molecule_id_t m_id,
    const float_t scheduled_time,
    const Reaction* unimol_rx
) {
  assert(unimol_rx != nullptr);
  assert(unimol_rx != (const Reaction*)-1);
  // the unimolecular reaction was already selected
  // FIXME: if there is more of them, mcell3 uses rng to select which to execute...
  Molecule& m = p.get_m(m_id);
  if (m.is_defunct()) {
    return;
  }

  // unimolecular reactions for surface molecules can be rescheduled,
  // ignore this action in this case
  if (scheduled_time != m.unimol_rx_time) {
    return;
  }

  assert(scheduled_time >= event_time && scheduled_time <= event_time + diffusion_time_step);
  int rx_res = outcome_unimolecular(p, m, scheduled_time - event_time, unimol_rx);
  assert(rx_res == RX_DESTROY);
}


// checks if reaction should probabilistically occur and if so,
// destroys reactants
// returns RX_DESTROY when reactants were destroyed
int DiffuseReactEvent::outcome_bimolecular(
    Partition& p,
    const Collision& collision,
    int path,
    float_t remaining_time_step
) {

  // might invalidate references

  int result =
      outcome_products_random(
        p,
        collision,
        remaining_time_step,
        path
      );

  if (result == RX_A_OK) {
    Molecule& reacA = p.get_m(collision.diffused_molecule_id);
    Molecule& reacB = p.get_m(collision.colliding_molecule_id);

#ifdef DEBUG_REACTIONS
    // reference printout first destroys B then A
    DUMP_CONDITION4(
      reacB.dump(p, "", "  defunct m:", world->get_current_iteration(), 0, false);
      reacA.dump(p, "", "  defunct m:", world->get_current_iteration(), 0, false);
    );
#endif

    // always for now
    // we used the reactants - remove them
    p.set_molecule_as_defunct(reacA);
    p.set_molecule_as_defunct(reacB);

    return RX_DESTROY;
  }


  return result;
}

// might return RX_BLOCKED if reaction cannot occur,
// returns 0 if positions were found
int DiffuseReactEvent::find_surf_product_positions(
    Partition& p,
    const Molecule* reacA,
    const Molecule* reacB,
    const Molecule* surf_reac,
    const Reaction* rx,
    small_vector<GridPos>& assigned_surf_product_positions) {

  uint num_surface_products = rx_util::get_num_surface_products(world, rx);

  small_vector<GridPos> recycled_surf_prod_positions; // this array contains information on where to place the surface products
  uint initiator_recycled_index = INDEX_INVALID;

  // find which tiles can be recycled
  if (reacA->is_surf()) {
    recycled_surf_prod_positions.push_back( GridPos::make_with_pos(p, *reacA) );
    if (reacA->id == surf_reac->id) {
      initiator_recycled_index = 0;
    }
  }

  if (reacB != nullptr && reacB->is_surf()) {
    recycled_surf_prod_positions.push_back( GridPos::make_with_pos(p, *reacB) );
    if (reacB->id == surf_reac->id) {
      initiator_recycled_index = recycled_surf_prod_positions.size() - 1;
    }
  }

  // do we need more tiles?
  TileNeighborVector vacant_neighbor_tiles;
  if (num_surface_products > recycled_surf_prod_positions.size()) {
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
    if (vacant_neighbor_tiles.size() + recycled_surf_prod_positions.size() < num_surface_products) {
      return RX_BLOCKED;
    }
  }

  // random assignment of positions
  uint num_tiles_to_recycle = min(rx->products.size(), recycled_surf_prod_positions.size());
  if (num_tiles_to_recycle == 1 && recycled_surf_prod_positions.size() >= 1) {
    // NOTE: this can be optimized - in case of a single product, it will just replace the initiator
    // and no initialization of recycled_surf_prod_positions is needed
    // there must be just one surface product for this variant
    assert(rx_util::get_num_surface_products(world, rx) == 1
        && "not sure, should not probably happen or this case is not handled");

    assert(initiator_recycled_index != INDEX_INVALID);
    assigned_surf_product_positions[0] = recycled_surf_prod_positions[initiator_recycled_index];
  }
  else if (num_tiles_to_recycle > 1) {
    uint next_available_index = 0;
    uint n_players = rx->get_num_players();
    uint n_reactants = rx->reactants.size();

    // assign recycled positions to products
    while (next_available_index < num_tiles_to_recycle) {
      // we must have the same number of random calls as in mcell3...
      uint rnd_num = rng_uint(&world->rng) % n_players;

      // continue until we got an index of a product
      if (rnd_num < n_reactants) {
        continue;
      }

      uint product_index = rnd_num - n_reactants;
      assert(product_index < rx->products.size());

      // we care only about surface molecules
      if (world->get_species(rx->products[product_index].species_id).is_vol()) {
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

    uint n_products = rx->products.size();
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
    float_t remaining_time_step,
    int path
) {
  assert(path == 0 && "Only single pathway is supported now");

  const Reaction* rx = collision.rx;
  assert(rx->reactants.size() == 1 || rx->reactants.size() == 2);

  Molecule* reacA = &p.get_m(collision.diffused_molecule_id);
  assert(reacA != nullptr);
  Molecule* reacB = nullptr;
  Molecule* surf_reac = nullptr;

  if (rx->reactants.size() == 2) {
    reacB = &p.get_m(collision.colliding_molecule_id);

    surf_reac = reacA->is_surf() ? reacA : reacB;

    /* Ensure that reacA and reacB are sorted in the same order as the rxn players. */
    /* Needed to maintain the same behavior as in mcell3 */
    if (SpeciesWithOrientation(reacA->species_id, reacA->get_orientation()) != rx->reactants[0]) {
      Molecule* tmp_mol = reacA;
      reacA = reacB;
      reacB = tmp_mol;
    }
    assert(rx->reactants[1].is_same_tolerate_orientation_none(reacB->species_id, reacB->get_orientation()));
  }
  else {
    surf_reac = reacA->is_surf() ? reacA : nullptr;
  }
  assert(rx->reactants[0].is_same_tolerate_orientation_none(reacA->species_id, reacA->get_orientation()));

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
    int res = find_surf_product_positions(p, reacA, reacB, surf_reac, rx, assigned_surf_product_positions);
    if (res == RX_BLOCKED) {
      return RX_BLOCKED;
    }
  }

  // free up tiles that we are probably going to reuse
  bool one_of_reactants_is_surf = false;
  if (reacA->is_surf()) {
    p.get_wall(reacA->s.wall_index).grid.reset_molecule_tile(reacA->s.grid_tile_index);
    reacA->s.grid_tile_index = TILE_INDEX_INVALID;
    one_of_reactants_is_surf = true;
  }

  if (reacB != nullptr && reacB->is_surf()) {
    p.get_wall(reacB->s.wall_index).grid.reset_molecule_tile(reacB->s.grid_tile_index);
    reacB->s.grid_tile_index = TILE_INDEX_INVALID;
    one_of_reactants_is_surf = true;
  }



  // create and place each product
  for (uint product_index = 0; product_index < rx->products.size(); product_index++) {
    const SpeciesWithOrientation& product = rx->products[product_index];
    const Species& species = world->get_species(product.species_id);

    molecule_id_t new_m_id;

    // set only for new vol mols when one of the reactants is surf, invalid by default
    WallTileIndexPair where_is_vm_created;

    float_t scheduled_time;
    if (rx->reactants.size() == 2 && species.is_vol() && !one_of_reactants_is_surf) {
      // bimolecular reaction
      // schedule new product for diffusion
      // collision.time is relative to the part that this molecule travels this diffusion step
      // so it needs to be scaled
      scheduled_time = event_time + diffusion_time_step - (remaining_time_step - collision.time * remaining_time_step);
    }
    else if (rx->reactants.size() == 2 && (species.is_surf() || one_of_reactants_is_surf)) {
      scheduled_time = event_time + collision.time;
    }
    else {
      // unimolecular reaction
      // reaction_time is the time when this new molecule was created
      scheduled_time = event_time + collision.time;
    }


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

        diffusion_util::tiny_diffuse_3D(p, new_vm, displacement, w.index, new_pos_after_diffuse);

        // update position and subpart if needed
        new_vm.v.pos = new_pos_after_diffuse;
        subpart_index_t new_subpart = p.get_subpartition_index(new_vm.v.pos);
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

      // get info on where to place the product
      const GridPos& new_grid_pos = assigned_surf_product_positions[product_index];
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

    // NOTE: in this time step, we will simply simulate all results of reactions regardless on the diffusion time step of the
    // particular product
    // we always create diffuse events, unimol react events are created elsewhere
    new_diffuse_or_unimol_react_actions.push_back(
        DiffuseOrUnimolReactionAction(
            DiffuseOrUnimolReactionAction::Type::DIFFUSE,
            new_m_id, scheduled_time,
            where_is_vm_created
    ));

  } // end for - product creation

  return RX_A_OK;
}

// ---------------------------------- unimolecular reactions ----------------------------------


int DiffuseReactEvent::outcome_unimolecular(
    Partition& p,
    Molecule& m,
    const float_t time_from_event_start,
    const Reaction* unimol_rx
) {
  molecule_id_t id = m.id;

  // creates new molecule(s) as output of the unimolecular reaction
  // !! might invalidate references (we might reorder defuncting and outcome call later)
  Collision collision(CollisionType::UNIMOLECULAR_VOLMOL, &p, m.id, time_from_event_start, m.v.pos, unimol_rx);
  //int outcome_res = outcome_products_random(p, unimol_rx, vm.v.pos, time_from_event_start, TIME_INVALID, 0);
  int outcome_res = outcome_products_random(p, collision, TIME_INVALID, 0);
  assert(outcome_res == RX_A_OK);

  // and defunct this molecule
  Molecule& m_new_ref = p.get_m(id);
#ifdef DEBUG_REACTIONS
  DUMP_CONDITION4(
    m_new_ref.dump(p, "", m_new_ref.is_vol() ? "Unimolecular vm defunct:" : "Unimolecular sm defunct:", world->get_current_iteration(), event_time + time_from_event_start, false);
  );
#endif
  p.set_molecule_as_defunct(m_new_ref);
  return RX_DESTROY;
}



// ---------------------------------- dumping methods ----------------------------------

void DiffuseReactEvent::dump(const std::string indent) {
  cout << indent << "Diffuse-react event:\n";
  std::string ind2 = indent + "  ";
  BaseEvent::dump(ind2);
  cout << ind2 << "diffusion_time_step: \t\t" << diffusion_time_step << " [float_t] \t\t\n";
}


} /* namespace mcell */
