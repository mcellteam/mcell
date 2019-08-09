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
#include <deque>
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
  COLLISION_INVALID,

  COLLISION_VOLMOL_VOLMOL,
  COLLISION_WALL_REDO,
  COLLISION_WALL_MISS,
  COLLISION_WALL_FRONT,
  COLLISION_WALL_BACK,

  COLLISION_SURFMOL_SURFMOL,

  COLLISION_VOLMOL_SURFMOL,

  COLLISION_UNIMOLECULAR_VOLMOL,
};


class collision_t;

#ifndef INDEXER_WA
typedef boost::container::small_vector<collision_t, 16> collision_vector_t;
typedef boost::container::flat_set<subpart_index_t> subpart_indices_set_t;
#else
typedef std::vector<collision_t> collision_vector_t;
typedef std::set<subpart_index_t> subpart_indices_set_t;
#endif

/**
 * Information about collision of 2 volume/surface molecules or a of a wall collision,
 * used in diffuse_react and in partition.
 */
class collision_t {
public:
  collision_t()
    : type(COLLISION_INVALID), partition(nullptr), diffused_molecule_id(MOLECULE_ID_INVALID), time(TIME_INVALID),
      colliding_molecule_id(MOLECULE_ID_INVALID), rx(nullptr), colliding_wall_index(WALL_INDEX_INVALID) {
  }

  // maybe create some static constructors with better names
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
    assert((type == COLLISION_VOLMOL_VOLMOL || type == COLLISION_VOLMOL_SURFMOL)
        && "This constructor must be used only for volmol or volsurf collisions");
  }

  collision_t(
      const collision_type_t type_,
      partition_t* partition_ptr,
      const molecule_id_t diffused_molecule_id_,
      const float_t time_, // time from event start
      const molecule_id_t colliding_molecule_id_,
      const reaction_t* rx_ptr
      )
    :
      type(type_),
      partition(partition_ptr),
      diffused_molecule_id(diffused_molecule_id_),
      time(time_),
      colliding_molecule_id(colliding_molecule_id_),
      rx(rx_ptr),
      colliding_wall_index(WALL_INDEX_INVALID) {
    assert(type == COLLISION_SURFMOL_SURFMOL && "This constructor must be used only for surfsurf collisions");
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
    assert((type == COLLISION_WALL_BACK || type == COLLISION_WALL_FRONT) && "This constructor must be used only for wall collisions");
  }

  collision_t(
      const collision_type_t type_,
      partition_t* partition_ptr,
      const molecule_id_t diffused_molecule_id_,
      const float_t time_,
      const vec3_t& pos_,
      const reaction_t* rx_ptr
      )
    :
      type(type_),
      partition(partition_ptr),
      diffused_molecule_id(diffused_molecule_id_),
      time(time_),
      pos(pos_),
      colliding_molecule_id(MOLECULE_ID_INVALID),
      rx(rx_ptr),
      colliding_wall_index(WALL_INDEX_INVALID) {
    assert(type == COLLISION_UNIMOLECULAR_VOLMOL && "This constructor must be used only for unimol volmol collisions");
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

  bool is_vol_mol_collision() const {
    return type == COLLISION_VOLMOL_VOLMOL;
  }

  bool is_wall_collision() const {
    assert(type != COLLISION_WALL_REDO && "Not sure yet what to do with redo");
    return type == COLLISION_WALL_FRONT || type == COLLISION_WALL_BACK;
  }

  void dump(partition_t& p, const std::string ind) const;
  std::string to_string() const;
  static void dump_array(partition_t& p, const collision_vector_t& vec);
};

/* contains information about the neighbors of the tile */
struct wall_tile_index_pair_t {
  wall_tile_index_pair_t(const wall_index_t wall_index_, const tile_index_t tile_index_)
    : wall_index(wall_index_), tile_index(tile_index_)
    {
  }

  wall_index_t wall_index;  /* surface grid the tile is on */
  tile_index_t tile_index;  /* index on that tile */
  //short int flag;           /* flag */ - not needed so far
};

struct tile_neighbor_vector_t: public std::deque<wall_tile_index_pair_t> {
  void dump(const std::string extra_comment, const std::string ind) {
    std::cout << ind << extra_comment << "\n";
    for (uint i = 0; i < size(); i++) {
      std::cout << ind << i << ": " << at(i).tile_index << "\n";
    }
  }
};

typedef small_vector<wall_index_t> wall_indices_t;

struct grid_pos_t {
  grid_pos_t()
    : initialized(false),
      wall_index(WALL_INDEX_INVALID), tile_index(TILE_INDEX_INVALID),
      pos_is_set(false), pos(POS_INVALID) {
  }

  static grid_pos_t make_with_pos(const partition_t& p, const molecule_t& sm) {
    assert(sm.is_surf());
    assert(sm.s.pos.is_valid());
    grid_pos_t res;

    res.initialized = true;
    res.wall_index = sm.s.wall_index;
    res.tile_index = sm.s.grid_tile_index;
    res.pos = sm.s.pos;
    res.pos_is_set = true;
    return res;
  }

  static grid_pos_t make_without_pos(const partition_t& p, const wall_tile_index_pair_t& wall_tile_index_pair) {
    assert(wall_tile_index_pair.tile_index != TILE_INDEX_INVALID);
    grid_pos_t res;

    res.initialized = true;
    res.wall_index = wall_tile_index_pair.wall_index;
    res.tile_index = wall_tile_index_pair.tile_index;
    res.pos_is_set = false;
    return res;
  }

  bool initialized; // was this info initialized
  wall_index_t wall_index;  /* wall where the tile is on */
  tile_index_t tile_index;  /* index on that tile */
  bool pos_is_set;
  vec2_t pos;
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
      const float_t time_up_to_event_end
  );

  // ---------------------------------- volume molecules ----------------------------------

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
      subpart_index_t& new_subpartition_index
  );

  bool collide_and_react_with_vol_mol(
      partition_t& p,
      collision_t& collision,
      vec3_t& displacement,
      float_t remaining_time_step,
      float_t r_rate_factor
  );

  int collide_and_react_with_surf_mol(
      partition_t& p,
      collision_t& collision,
      float_t remaining_time_step,
      float_t r_rate_factor,
      float_t current_molecule_time
  );

  // ---------------------------------- surface molecules ----------------------------------

  void diffuse_surf_molecule(
      partition_t& p,
      const molecule_id_t sm_id,
      const float_t current_time,
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

  bool react_2D_all_neighbors(
      partition_t& p,
      molecule_t& sm,
      const float_t current_time,
      const float_t remaining_time_step
  );

  // ---------------------------------- reactions ----------------------------------

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
      collision_t& collision,
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
