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

#include <iostream>
#include <sstream>
#include <algorithm>
#include <boost/container/flat_set.hpp>

extern "C" {
#include "rng.h"
#include "mcell_structs.h"
#include "logging.h"
}

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

using namespace std;

namespace mcell {

void diffuse_react_event_t::step() {
  assert(world->partitions.size() == 1 && "Must extend cache to handle multiple partitions");

  // for each partition
  for (partition_t& p: world->partitions) {
    // diffuse molecules from volume_molecule_indices_per_time_step that have the current diffusion_time_step
    uint32_t time_step_index = p.get_molecule_list_index_for_time_step(diffusion_time_step);
    if (time_step_index != TIME_STEP_INDEX_INVALID) {
      diffuse_molecules(p, p.get_volume_molecule_ids_for_time_step_index(time_step_index));
    }
  }
}


void diffuse_react_event_t::diffuse_molecules(partition_t& p, const std::vector<molecule_id_t>& molecule_ids) {
  float_t event_time_end = event_time + diffusion_time_step;

  // we need to stricly follow the ordering in mcell3, therefore steps 2) and 3) do not use the time
  // for which they were scheduled but rather simply the order in which these "microevents" were created

  // 1) first diffuse already existing molecules
  uint32_t existing_mols_count = molecule_ids.size();
  for (uint32_t i = 0; i < existing_mols_count; i++) {
    molecule_id_t id = molecule_ids[i];
    // existing molecules - simulate whole time step
    diffuse_single_molecule(p, id, diffusion_time_step, event_time_end);
  }

  // 2) we need to take care of unimolecular reactions that were scheduled for this time step
  // in the previous time steps; these unimolecular reaction microevents are handled like they are in a queue

  // first get the calndar owned by partition that contains actions that are scheduled
  uint32_t time_step_index = p.get_or_add_molecule_list_index_for_time_step(diffusion_time_step);
  partition_t::calendar_for_unimol_rxs_t& calendar_for_unimol_rxs =
      p.get_unimolecular_actions_calendar_for_time_step_index(time_step_index);

  // then get bucket for this this->time corresponding to our event time
  uint64_t bucket_index = calendar_for_unimol_rxs.get_bucket_index_for_time(event_time);
  if (bucket_index != BUCKET_INDEX_INVALID) {

    // FIXME: somehow simplify this
    fifo_bucket_t<diffuse_or_unimol_react_action_t>& bucket = calendar_for_unimol_rxs.get_bucket_with_index(bucket_index);
    std::vector<diffuse_or_unimol_react_action_t>& actions = bucket.events;
    for (const diffuse_or_unimol_react_action_t& unimol_action: actions) {
      react_unimol_single_molecule(p, unimol_action.id, unimol_action.scheduled_time, unimol_action.unimol_rx);
    }

    // remove bucket and also all the older ones from out internal scheduler because we processed it
    // FIXME: we can remove multiple items at once
    uint64_t i = 0;
    while (i <= bucket_index) {
      calendar_for_unimol_rxs.pop_bucket();
      i++;
    }
  }

  // 3) simulate remaining time of molecules created with reactions
  // need to call .size() each iteration because the size can increase,
  // again, we are using it as a queue and we do not follow the time when they were created
  for (uint32_t i = 0; i < new_diffuse_or_unimol_react_actions.size(); i++) {
    const diffuse_or_unimol_react_action_t& action = new_diffuse_or_unimol_react_actions[i];

    if (action.type == diffuse_or_unimol_react_action_t::DIFFUSE) {
      diffuse_single_molecule(p, action.id, diffusion_time_step - action.scheduled_time, event_time_end);
    }
    else {
      react_unimol_single_molecule(p, action.id, action.scheduled_time, action.unimol_rx);
    }
  }

  new_diffuse_or_unimol_react_actions.clear();
}


// get displacement based on scale (related to diffusion constant) and gauss random number
static void pick_displacement(float_t scale, rng_state& rng, vec3_t& displacement) {
  displacement.x = scale * rng_gauss(&rng) * .70710678118654752440;
  displacement.y = scale * rng_gauss(&rng) * .70710678118654752440;
  displacement.z = scale * rng_gauss(&rng) * .70710678118654752440;
}


// determine how far will our diffused molecule move
static void compute_displacement(
    species_t& sp,
    rng_state& rng,
    float_t remaining_time_step,
    vec3_t& displacement,
    float_t& r_rate_factor) {

  float_t rate_factor = (remaining_time_step == 1.0) ? 1.0 : sqrt(remaining_time_step);
  r_rate_factor = 1.0 / rate_factor;
  pick_displacement(sp.space_step * rate_factor, rng, displacement);
}


void diffuse_react_event_t::diffuse_single_molecule(
    partition_t& p,
    const molecule_id_t vm_id,
    const float_t time_up_to_event_end,
    const float_t event_time_end
) {

  volume_molecule_t& vm = p.get_vm(vm_id);

  if (vm.is_defunct())
    return;

  // if the molecule is a "newbie", its unimolecular reaction was not yet scheduled
  if ((vm.flags & ACT_NEWBIE) != 0) {
    create_unimol_rx_action(p, vm, time_up_to_event_end);
    vm.flags &= ~ACT_NEWBIE;
  }


#ifdef DEBUG_DIFFUSION
  DUMP_CONDITION4(
    // the subtraction of diffusion_time_step doesn't make much sense but is needed to make the dump the same as in mcell3
    // need to check it further
    vm.dump(world, "", "Diffusing vm:", world->current_iteration, event_time_end - time_up_to_event_end - diffusion_time_step);
  );
#endif

  // we might need to adjust remaining time step if this molecule has a unimolecular reaction
  // within this event's time step range
  float_t remaining_time_step;
  if (vm.unimol_rx_time < event_time_end) { // unlikely
    assert(vm.unimol_rx_time >= event_time && "Missed unimol rx");

    // the value of remaining_time_step passed as argument is time up to event_time_end
    float_t prev_time_from_event_start = diffusion_time_step - time_up_to_event_end;
    float_t new_time_from_event_start = vm.unimol_rx_time - event_time;
    assert(new_time_from_event_start >= prev_time_from_event_start && "Unimol rx cannot be scheduled to the past");

    remaining_time_step = new_time_from_event_start - prev_time_from_event_start;
  }
  else {
    remaining_time_step = time_up_to_event_end;
  }

  species_t& species = world->species[vm.species_id];

  // diffuse each molecule - get information on position change
  // TBD: reflections
  vec3_t displacement;
  float_t r_rate_factor;
  compute_displacement(species, world->rng, remaining_time_step, displacement, r_rate_factor);

#ifdef DEBUG_DIFFUSION
  DUMP_CONDITION4(
    displacement.dump("  displacement:", "");
  );
#endif
  // note: we are ignoring use_expanded_list setting compared to mcell3

  // detect collisions with other molecules
  vec3_t remaining_displacement = displacement;
  uint32_t new_subpart_index;
  vec3_t new_pos;
  ray_trace_state_t state;
  collision_vector_t molecule_collisions;
  bool was_defunct = false;
  wall_index_t reflected_wall_index = WALL_INDEX_INVALID;
  do {
    state =
        ray_trace(
            p, vm /* changes position */,
            reflected_wall_index,
            remaining_displacement,
            molecule_collisions,
            new_pos,
            new_subpart_index
        );

    // sort current collisions by time
    sort( molecule_collisions.begin(), molecule_collisions.end(),
        [ ]( const auto& lhs, const auto& rhs )
        {
          if (lhs.time < rhs.time) {
            return true;
          }
          else if (lhs.time > rhs.time) {
            return false;
          }
          else if (lhs.type == COLLISION_VOLMOL_VOLMOL && rhs.type == COLLISION_VOLMOL_VOLMOL) {
            // mcell3 returns collisions with molecules ordered descending by the molecule index
            // we need to maintain this behavior (needed only for identical results)
            return lhs.colliding_molecule_id > rhs.colliding_molecule_id;
          }
          else {
            return false;
          }
        }
    );

  #ifdef DEBUG_COLLISIONS
    DUMP_CONDITION4(
      collision_t::dump_array(p, molecule_collisions);
    );
  #endif

    // evaluate and possible execute collisions and reactions
    for (size_t collision_index = 0; collision_index < molecule_collisions.size(); collision_index++) {
      collision_t& collision = molecule_collisions[collision_index];

      assert(collision.time >= 0 && collision.time <= 1);

      if (collision.type == COLLISION_VOLMOL_VOLMOL) {
        // ignoring immediate collisions
        if (collision.time < EPS) {
            continue;
        }

        // evaluate reaction associated with this collision
        // for now, do the change right away, but we might need to cache these changes and
        // do them after all diffusions were finnished
        // warning: might invalidate references to p.volume_molecules array! returns true in that case
        if (collide_and_react_with_vol_mol(p, collision, remaining_displacement, remaining_time_step, r_rate_factor)) {
          // molecule was destroyed
           was_defunct = true;
          break;
        }
      }
      else if (collision.is_wall_collision()) {
        volume_molecule_t& vm_new_ref = p.get_vm(vm_id);
        int res = coll_util::reflect_or_periodic_bc(p, collision, vm_new_ref, remaining_displacement, remaining_time_step, reflected_wall_index);
        assert(res == 0 && "Periodic box BCs are not supported yet");

        // molecule could have been moved
        subpart_index_t subpart_after_wall_hit = p.get_subpartition_index(vm_new_ref.pos);
        // change subpartition if needed
        p.change_molecule_subpartition(vm_new_ref, subpart_after_wall_hit);


        break; // we reflected, do ray_trace again
      }
    }

  } while (unlikely(state != RAY_TRACE_FINISHED && !was_defunct));


  if (!was_defunct) {
    // need to get a new reference
    volume_molecule_t& vm_new_ref = p.get_vm(vm_id);

    // finally move molecule to its destination
    vm_new_ref.pos = new_pos;

    // are we still in the same partition or do we need to move?
    bool move_to_another_partition = !p.in_this_partition(vm_new_ref.pos);
    if (move_to_another_partition) {
      mcell_log("Error: Crossing partitions is not supported yet.\n");
      exit(1);
    }

    // change subpartition
    p.change_molecule_subpartition(vm_new_ref, new_subpart_index);
  }
}


// collect possible collisions for molecule vm that has to displace by remaining_displacement,
// returns possible collisions in molecule_collisions, new position in new_pos and
// index of the new subparition in new_subpart_index
// later, this will check collisions until a wall is hit
// we assume that wall collisions do not occur so often
ray_trace_state_t diffuse_react_event_t::ray_trace(
    partition_t& p,
    volume_molecule_t& vm, // molecule that we are diffusing, we are changing its pos  and possibly also subvolume
    const wall_index_t previous_reflected_wall, // is WALL_INDEX_INVALID when our molecule did not replect from anything this iddfusion step yet
    vec3_t& remaining_displacement, // in/out - recomputed if there was a reflection
    collision_vector_t& collisions, // both mol mol and wall collisions
    vec3_t& new_pos,
    uint32_t& new_subpart_index
    ) {
  ray_trace_state_t res_state = RAY_TRACE_FINISHED;
  collisions.clear();

  // first get what subpartitions might be relevant
  subpart_indices_vector_t crossed_subparts_for_walls;
  subpart_indices_set_t crossed_subparts_for_molecules;
  uint32_t last_subpartition_index;
  coll_util::collect_crossed_subparts(
      p, vm, remaining_displacement,
      world->world_constants.rx_radius_3d,  world->world_constants.subpartition_edge_length,
      true, crossed_subparts_for_walls,
      crossed_subparts_for_molecules, last_subpartition_index
  );


  vec3_t& displacement_up_to_wall_collision = remaining_displacement;

  // check wall collisions in the crossed subparitions,
  // stop at first crossing because crossed_subparts_for_walls are ordered
  // and we are sure that if we hit a wall in the actual supartition, we cannot
  // possibly hit another wall in a subparition that follows
  if (!crossed_subparts_for_walls.empty()) {
    for (uint32_t subpart_index: crossed_subparts_for_walls) {

      coll_util::collect_wall_collisions( // mcell3 does this onlty for the current subvol
          p,
          vm,
          subpart_index,
          previous_reflected_wall,
          world->rng,
          displacement_up_to_wall_collision,
          collisions
      );

      // hits are simply mixed with reactions... -
      // TODO: use the same structure, also might cache walls that can collide
  /*
      if (!collisions.empty()) {
        res_state = RAY_TRACE_HIT_WALL;

        // TODO: optimization?
        // displacement_up_to_wall_collision -

        // we do not care about wall collisions in following subvolumes
        break;
      }*/
    }
    if (!collisions.empty()) {
      res_state = RAY_TRACE_HIT_WALL;
    }
  }

  if (res_state == RAY_TRACE_HIT_WALL) {
    // recompute collect_crossed_subparts if there was a wall collision
    crossed_subparts_for_molecules.clear();
    coll_util::collect_crossed_subparts(
        p, vm, displacement_up_to_wall_collision,
        world->world_constants.rx_radius_3d,
        world->world_constants.subpartition_edge_length,
        false, crossed_subparts_for_walls, // not filled this time
        crossed_subparts_for_molecules, last_subpartition_index
    );
  }

  float_t radius = world->world_constants.rx_radius_3d;

  // check molecule collisions for each SP
  for (uint32_t subpart_index: crossed_subparts_for_molecules) {
    // get cached reacting molecules for this SP
    subpartition_mask_t& sp_reactants = p.get_volume_molecule_reactants(subpart_index, vm.species_id);

    // for each molecule in this SP
    for (uint32_t colliding_vm_index: sp_reactants) {
      coll_util::collide_mol_loop_body(
          world,
          p,
          vm,
          colliding_vm_index,
          displacement_up_to_wall_collision,
          radius,
          collisions
      );
    }
  }

  // the value is valid only when RAY_TRACE_FINISHED is returned
  new_subpart_index = last_subpartition_index; // FIXME: this is valid only when there was no collision
  new_pos = vm.pos + remaining_displacement; // FIXME: this i overwritten in case of a collision, should I do somthing about it?

  return res_state; // no wall was hit
}


// handle collision of two volume molecules: checks probability of reaction,
// executes this reaction, removes reactants and creates products
// returns true if reaction has occured and the first reactant was destroyed
bool diffuse_react_event_t::collide_and_react_with_vol_mol(
    partition_t& p,
    collision_t& collision,
    vec3_t& displacement,
    float_t remaining_time_step,
    float_t r_rate_factor
)  {

  volume_molecule_t& colliding_molecule = p.get_vm(collision.colliding_molecule_id); // am
  volume_molecule_t& diffused_molecule = p.get_vm(collision.diffused_molecule_id); // m

  // returns 1 when there are no walls at all
  //TBD: double factor = exact_disk(
  float_t factor = exact_disk_util::exact_disk(
      p, collision.pos, displacement, p.get_world_constants().rx_radius_3d,
      diffused_molecule, colliding_molecule,
      p.get_world_constants().use_expanded_list
  );

  if (factor < 0) { /* Probably hit a wall, might have run out of memory */
    return 0; /* Reaction blocked by a wall */
  }

  const reaction_t& rx = *collision.rx;
  //  rx->prob_t is always NULL in out case update_probs(world, rx, m->t);
  // returns which reaction pathway to take
  float_t scaling = factor * r_rate_factor;
  int i = rx_util::test_bimolecular(
    rx, world->rng, colliding_molecule, diffused_molecule, scaling);

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


void diffuse_react_event_t::create_unimol_rx_action(
    partition_t& p,
    volume_molecule_t& vm,
    float_t remaining_time_step
) {
  float_t curr_time = event_time + diffusion_time_step - remaining_time_step;
  assert(curr_time >= 0);

  const reaction_t* rx = rx_util::pick_unimol_rx(world, vm.species_id);
  if (rx == nullptr) {
    return;
  }

  float_t time_from_now = rx_util::compute_unimol_lifetime(world, vm, rx);

  float_t scheduled_time = curr_time + time_from_now;

  // we need to store the end time to the molecule because oit is needed in diffusion to
  // figure out whether we should do the whole time step
  vm.unimol_rx_time = scheduled_time;

  // now, there are two queues - local for this timestep
  // and global in partition for the following timesteps7
  diffuse_or_unimol_react_action_t unimol_react_action(
      vm.id, scheduled_time, diffuse_or_unimol_react_action_t::UNIMOL_REACT, rx);

  if (scheduled_time < event_time + diffusion_time_step) {
    // handle this iteration
    new_diffuse_or_unimol_react_actions.push_back(unimol_react_action);
  }
  else {
    p.add_unimolecular_action(diffusion_time_step, unimol_react_action);
  }
}


// based on mcell3's check_for_unimolecular_reaction
// might invalidate vm references
void diffuse_react_event_t::react_unimol_single_molecule(
    partition_t& p,
    const molecule_id_t vm_id,
    const float_t scheduled_time,
    const reaction_t* unimol_rx
) {
  assert(unimol_rx != nullptr);
  // the unimolecular reaction was already selected
  // FIXME: if there is more of them, mcell3 uses rng to select which to execute...
  volume_molecule_t& vm = p.get_vm(vm_id);
  if (vm.is_defunct()) {
    return;
  }

  assert(scheduled_time >= event_time && scheduled_time <= event_time + diffusion_time_step);
  int rx_res = outcome_unimolecular(p, vm, scheduled_time - event_time, unimol_rx);
  assert(rx_res == RX_DESTROY);
}


// checks if reaction should probabilistically occur and if so,
// destroys reactants
// returns RX_DESTROY when reactants were destroyed
int diffuse_react_event_t::outcome_bimolecular(
    partition_t& p,
    collision_t& collision,
    int path,
    float_t remaining_time_step
) {

  // might invalidate references!
  int result =
      outcome_products_random(
        p,
        collision.rx,
        collision.pos,
        collision.time,
        remaining_time_step,
        path
      );

  if (result == RX_A_OK) {
    volume_molecule_t& reacA = p.get_vm(collision.diffused_molecule_id);
    volume_molecule_t& reacB = p.get_vm(collision.colliding_molecule_id);

#ifdef DEBUG_REACTIONS
    // reference printout first destroys B then A
    DUMP_CONDITION4(
      reacB.dump(world, "", "  defunct vm:", world->current_iteration);
      reacA.dump(world, "", "  defunct vm:", world->current_iteration);
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


// why is this called "random"? - check if reaction occurs is in test_bimolecular
// mcell3 version returns  cross_wall ? RX_FLIP : RX_A_OK;
// ! might invalidate references
int diffuse_react_event_t::outcome_products_random(
    partition_t& p,
    const reaction_t* rx,
    const vec3_t& pos,
    float_t reaction_time,
    float_t remaining_time_step,
    int path
) {
  assert(path == 0 && "Only single pathway is supported now");
  // we can have just one product for now and no walls

  // create and place each product
  for (const species_with_orientation_t& product: rx->products) {
    volume_molecule_t vm(MOLECULE_ID_INVALID, product.species_id, pos);

    volume_molecule_t& new_vm = p.add_volume_molecule(vm, world->species[vm.species_id].time_step);
    new_vm.flags =  ACT_NEWBIE | TYPE_VOL | IN_VOLUME | ACT_DIFFUSE;

  #ifdef DEBUG_REACTIONS
    DUMP_CONDITION4(
      new_vm.dump(world, "", "  created vm:", world->current_iteration);
    );
  #endif

    float_t scheduled_time;
    if (rx->reactants.size() == 2) {
      // bimolecular reaction
      // schedule new product for diffusion
      // collision.time is relative to the part that this molecule travels this diffusion step
      // so it needs to be scaled
      scheduled_time = diffusion_time_step - (remaining_time_step - reaction_time*remaining_time_step);
    }
    else {
      // unimolecular reaction
      // reaction_time is the time when this new molecule was created
      scheduled_time = reaction_time;
    }

    // NOTE: in this time step, we will simply simulate all results of reactions regardless on the diffusion time step of the
    // particular product
    // we alway create diffuse events, unimol react events are created elsewhere
    new_diffuse_or_unimol_react_actions.push_back(
        diffuse_or_unimol_react_action_t(new_vm.id, scheduled_time, diffuse_or_unimol_react_action_t::DIFFUSE));
  }
  return RX_A_OK;
}

// ---------------------------------- unimolecular reactions ----------------------------------


int diffuse_react_event_t::outcome_unimolecular(
    partition_t& p,
    volume_molecule_t& vm,
    const float_t time_from_event_start,
    const reaction_t* unimol_rx
) {
  molecule_id_t id = vm.id;

  // creates new molecule(s) as output of the unimolecular reaction
  // !! might invalidate references (we might reorder defuncting and outcome call later)
  int outcome_res = outcome_products_random(p, unimol_rx, vm.pos, time_from_event_start, TIME_INVALID, 0);
  assert(outcome_res == RX_A_OK);

  // and defunct this molecule
  volume_molecule_t& vm_new_ref = p.get_vm(id);
#ifdef DEBUG_REACTIONS
  DUMP_CONDITION4(
    vm_new_ref.dump(world, "", "Unimolecular vm defunct:", world->current_iteration, time_from_event_start);
  );
#endif
  p.set_molecule_as_defunct(vm_new_ref);
  return RX_DESTROY;
}



// ---------------------------------- dumping methods ----------------------------------

void diffuse_react_event_t::dump(const std::string indent) {
  cout << indent << "Diffuse-react event:\n";
  std::string ind2 = indent + "  ";
  base_event_t::dump(ind2);
  cout << ind2 << "diffusion_time_step: \t\t" << diffusion_time_step << " [float_t] \t\t\n";
}


void collision_t::dump(partition_t& p, const std::string ind) const {
  cout << ind << "diffused_molecule:\n";
  p.get_vm(diffused_molecule_id).dump(ind + "  ");
  if (type == COLLISION_VOLMOL_VOLMOL) {
    cout << ind << "colliding_molecule:\n";
    p.get_vm(colliding_molecule_id).dump(ind + "  ");
    cout << ind << "reaction:";
    rx->dump(ind + "  ");
  }
  else {
    cout << ind << "colliding_wall_index: " << colliding_wall_index << "\n";
  }

  cout << "time: \t\t" << time << " [float_t] \t\t\n";
  cout << "position: \t\t" << pos << " [vec3_t] \t\t\n";
}


string collision_t::to_string() const {
  stringstream ss;
  if (type == COLLISION_VOLMOL_VOLMOL) {
    ss << "coll_idx: " << colliding_molecule_id;
  }
  else {
    //ss << "wall_idx: " << colliding_wall_index;
  }

  ss << ", time: " << time << ", pos: " << pos;
  return ss.str();
}


void collision_t::dump_array(partition_t& p, const collision_vector_t& vec) {
  // printed in resverse - same as
  for (size_t i = 0; i < vec.size(); i++) {
    const char* str_type = (vec[i].type == COLLISION_VOLMOL_VOLMOL) ? "mol collision " : "wall collision ";
    cout << "  " << str_type << i << ": " << vec[i].to_string() << "\n";
  }
}

} /* namespace mcell */
