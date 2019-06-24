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

#ifndef SRC4_DIFFUSE_REACT_EVENT_H_
#define SRC4_DIFFUSE_REACT_EVENT_H_

#include <vector>
#include "../libs/boost/container/small_vector.hpp"

#include "base_event.h"
#include "partition.h"
#include "reaction.h"


namespace mcell {

class partition_t;
class molecule_t;
class species_t;


enum ray_trace_state_t {
  RAY_TRACE_HIT_UNDEFINED,
  RAY_TRACE_HIT_SUBPARTITION,
  RAY_TRACE_HIT_WALL,
  RAY_TRACE_FINISHED
};


enum collision_type_t {
  COLLISION_VOLMOL_VOLMOL,
  COLLISION_WALL_REDO,
  COLLISION_WALL_MISS,
  COLLISION_WALL_FRONT,
  COLLISION_WALL_BACK
};


class collision_t;
typedef boost::container::small_vector<collision_t, 16> collision_vector_t;
typedef boost::container::flat_set<subpart_index_t> subpart_indices_set_t;


/**
 * Information about collision of 2 volume molecules or a of a wall collision,
 * used in diffuse_react and in partition.
 */
class collision_t {
public:
  collision_t(
      const collision_type_t type_,
      partition_t* partition_ptr,
      const molecule_id_t diffused_molecule_id_,
      const float_t time_,
      const vec3_t& pos_,
      const molecule_id_t colliding_molecule_id_,
      const reaction_t* rx_ptr
      )
    :
      type(type_),
      partition(partition_ptr),
      diffused_molecule_id(diffused_molecule_id_),
      time(time_),
      pos(pos_),
      colliding_molecule_id(colliding_molecule_id_),
      rx(rx_ptr),
      colliding_wall_index(WALL_INDEX_INVALID) {
    assert(type == COLLISION_VOLMOL_VOLMOL && "This constructor must be used only for volmol collisions");
  }

  collision_t(
      const collision_type_t type_,
      partition_t* partition_ptr,
      const molecule_id_t diffused_molecule_id_,
      const float_t time_,
      const vec3_t& pos_,
      const wall_index_t colliding_wall_index_
      )
    :
      type(type_),
      partition(partition_ptr),
      diffused_molecule_id(diffused_molecule_id_),
      time(time_),
      pos(pos_),
      colliding_molecule_id(MOLECULE_ID_INVALID),
      rx(nullptr),
      colliding_wall_index(colliding_wall_index_) {
    assert(type != COLLISION_VOLMOL_VOLMOL && "This constructor must be used only for wall collisions");
  }

  collision_type_t type;
  partition_t* partition;
  molecule_id_t diffused_molecule_id;
  float_t time;
  vec3_t pos;

  // valid only for COLLISION_VOLMOL_VOLMOL
  molecule_id_t colliding_molecule_id;
  const reaction_t* rx;

  // valid only for COLLISION_WALL*
  wall_index_t colliding_wall_index;

  bool is_wall_collision() const {
    assert(type != COLLISION_WALL_REDO && "Not sure yet what to do with redo");
    return type == COLLISION_WALL_FRONT || type == COLLISION_WALL_BACK;
  }

  void dump(partition_t& p, const std::string ind) const;
  std::string to_string() const;
  static void dump_array(partition_t& p, const collision_vector_t& vec);
};



/**
 * Diffuse all molecules with a given time step.
 * When a molecule is diffused, it is checked for collisions and reactions
 * are evaluated. Molecules newly created in reactions and diffused with their
 * remaining time.
 */
class diffuse_react_event_t : public base_event_t {
public:
  diffuse_react_event_t(world_t* world_, const float_t diffusion_time_step_) :
    base_event_t(EVENT_TYPE_INDEX_DIFFUSE_REACT),
    world(world_),
    diffusion_time_step(diffusion_time_step_) {

    // repeat this event each time step
    periodicity_interval = diffusion_time_step;
  }
  void step();
  void dump(const std::string indent);

  world_t* world;

  // this event diffuses all molecules that have this diffusion time_step
  float_t diffusion_time_step;

private:
  // molecules newly created in reactions
  std::vector<diffuse_or_unimol_react_action_t> new_diffuse_or_unimol_react_actions;

  void diffuse_molecules(partition_t& p, const std::vector< molecule_id_t >& indices);

  void diffuse_single_molecule(
      partition_t& p,
      const molecule_id_t vm_id,
      const float_t time_up_to_event_end,
      const float_t event_time_end
  );

  void diffuse_vol_molecule(
      partition_t& p,
      const molecule_id_t vm_id,
      const float_t remaining_time_step,
      bool& was_defunct,
      vec3_t& new_pos,
      subpart_index_t& new_subpart_index
  );

  ray_trace_state_t ray_trace_vol(
      partition_t& p,
      molecule_t& vm, // molecule that we are diffusing, we are changing its pos  and possibly also subvolume
      const wall_index_t previous_reflected_wall, // is WALL_INDEX_INVALID when our molecule did not replect from anything this iddfusion step yet
      vec3_t& remaining_displacement, // in/out - recomputed if there was a reflection
      collision_vector_t& molecule_collisions, // possible reactions in this part of way marching, ordered by time
      vec3_t& new_position,
      uint32_t& new_subpartition_index
  );

  bool collide_and_react_with_vol_mol(
      partition_t& p,
      collision_t& collision,
      vec3_t& displacement,
      float_t remaining_time_step,
      float_t r_rate_factor
  );

  void diffuse_surf_molecule(
      partition_t& p,
      const molecule_id_t sm_id,
      const float_t remaining_time_step,
      bool& was_defunct,
      vec2_t& new_loc,
      wall_index_t& new_wall_index,
      float_t& advance_time
  );

  wall_index_t ray_trace_surf(
      partition_t& p,
      const species_t& species,
      molecule_t& sm,
      vec2_t& remaining_displacement,
      vec2_t& new_pos,
      bool& was_defunct
  );

  int outcome_bimolecular(
      partition_t& p,
      collision_t& collision,
      int path,
      float_t remaining_time_step
  );

  int outcome_unimolecular(
      partition_t& p,
      molecule_t& vm,
      const float_t scheduled_time,
      const reaction_t* unimol_rx
  );

  int outcome_products_random(
      partition_t& p,
      const reaction_t* rx,
      const vec3_t& pos,
      float_t reaction_time,
      float_t remaining_time_step,
      int path
  );

  void create_unimol_rx_action(
      partition_t& p,
      molecule_t& vm,
      float_t remaining_time_step
  );

  void react_unimol_single_molecule(
      partition_t& p,
      const molecule_id_t vm_id,
      const float_t scheduled_time,
      const reaction_t* unimol_rx
  );
};

} // namespace mcell

#endif // SRC4_DIFFUSE_REACT_EVENT_H_
