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


namespace MCell {

class World;
class Partition;
class Molecule;
class Species;


enum class RayTraceState {
  UNDEFINED,
  HIT_SUBPARTITION,
  RAY_TRACE_HIT_WALL,
  FINISHED
};


enum class CollisionType {
  INVALID,

  VOLMOL_VOLMOL,
  WALL_REDO,
  WALL_MISS,
  WALL_FRONT,
  WALL_BACK,

  SURFMOL_SURFMOL,

  VOLMOL_SURFMOL,

  UNIMOLECULAR_VOLMOL,
};


class Collision;

typedef boost::container::small_vector<Collision, 16> collision_vector_t;
typedef boost::container::flat_set<subpart_index_t> subpart_indices_set_t;

/**
 * Information about collision of 2 volume/surface molecules or a of a wall collision,
 * used in diffuse_react and in partition.
 */
// TODO: move to collision utils?
class Collision {
public:
  Collision()
    : type(CollisionType::INVALID), partition(nullptr), diffused_molecule_id(MOLECULE_ID_INVALID), time(TIME_INVALID),
      colliding_molecule_id(MOLECULE_ID_INVALID), rx(nullptr), colliding_wall_index(WALL_INDEX_INVALID) {
  }

  // maybe create some static constructors with better names
  Collision(
      const CollisionType type_,
      Partition* partition_ptr,
      const molecule_id_t diffused_molecule_id_,
      const float_t time_,
      const vec3_t& pos_,
      const molecule_id_t colliding_molecule_id_,
      const Reaction* rx_ptr
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
    assert((type == CollisionType::VOLMOL_VOLMOL || type == CollisionType::VOLMOL_SURFMOL)
        && "This constructor must be used only for volmol or volsurf collisions");
  }

  Collision(
      const CollisionType type_,
      Partition* partition_ptr,
      const molecule_id_t diffused_molecule_id_,
      const float_t time_, // time from event start
      const molecule_id_t colliding_molecule_id_,
      const Reaction* rx_ptr
      )
    :
      type(type_),
      partition(partition_ptr),
      diffused_molecule_id(diffused_molecule_id_),
      time(time_),
      colliding_molecule_id(colliding_molecule_id_),
      rx(rx_ptr),
      colliding_wall_index(WALL_INDEX_INVALID) {
    assert(type == CollisionType::SURFMOL_SURFMOL && "This constructor must be used only for surfsurf collisions");
  }

  Collision(
      const CollisionType type_,
      Partition* partition_ptr,
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
    assert((type == CollisionType::WALL_BACK || type == CollisionType::WALL_FRONT) && "This constructor must be used only for wall collisions");
  }

  Collision(
      const CollisionType type_,
      Partition* partition_ptr,
      const molecule_id_t diffused_molecule_id_,
      const float_t time_,
      const vec3_t& pos_,
      const Reaction* rx_ptr
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
    assert(type == CollisionType::UNIMOLECULAR_VOLMOL && "This constructor must be used only for unimol volmol collisions");
  }


  CollisionType type;
  Partition* partition;
  molecule_id_t diffused_molecule_id;
  float_t time;
  vec3_t pos;

  // valid only for COLLISION_VOLMOL_VOLMOL
  molecule_id_t colliding_molecule_id;
  const Reaction* rx;

  // valid only for COLLISION_WALL*
  wall_index_t colliding_wall_index;

  bool is_vol_mol_collision() const {
    return type == CollisionType::VOLMOL_VOLMOL;
  }

  bool is_wall_collision() const {
    assert(type != CollisionType::WALL_REDO && "Not sure yet what to do with redo");
    return type == CollisionType::WALL_FRONT || type == CollisionType::WALL_BACK;
  }

  void dump(Partition& p, const std::string ind) const;
  std::string to_string(const Partition& p) const;
  static void dump_array(Partition& p, const collision_vector_t& vec);
};

/* contains information about the neighbors of the tile */
class WallTileIndexPair {
public:
  WallTileIndexPair(const wall_index_t wall_index_, const tile_index_t tile_index_)
    : wall_index(wall_index_), tile_index(tile_index_)
    {
  }

  wall_index_t wall_index;  /* surface grid the tile is on */
  tile_index_t tile_index;  /* index on that tile */
  //short int flag;           /* flag */ - not needed so far
};

class TileNeighborVector: public std::deque<WallTileIndexPair> {
public:
  void dump(const std::string extra_comment, const std::string ind) {
    std::cout << ind << extra_comment << "\n";
    for (uint i = 0; i < size(); i++) {
      std::cout << ind << i << ": " << at(i).tile_index << "\n";
    }
  }
};

typedef small_vector<wall_index_t> wall_indices_t;

class GridPos {
public:
  GridPos()
    : initialized(false),
      wall_index(WALL_INDEX_INVALID), tile_index(TILE_INDEX_INVALID),
      pos_is_set(false), pos(POS_INVALID) {
  }

  static GridPos make_with_pos(const Partition& p, const Molecule& sm) {
    assert(sm.is_surf());
    assert(sm.s.pos.is_valid());
    GridPos res;

    res.initialized = true;
    res.wall_index = sm.s.wall_index;
    res.tile_index = sm.s.grid_tile_index;
    res.pos = sm.s.pos;
    res.pos_is_set = true;
    return res;
  }

  static GridPos make_without_pos(const Partition& p, const WallTileIndexPair& wall_tile_index_pair) {
    assert(wall_tile_index_pair.tile_index != TILE_INDEX_INVALID);
    GridPos res;

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
class DiffuseReactEvent : public BaseEvent {
public:
  DiffuseReactEvent(World* world_, const float_t diffusion_time_step_) :
    BaseEvent(EVENT_TYPE_INDEX_DIFFUSE_REACT),
    world(world_),
    diffusion_time_step(diffusion_time_step_) {

    // repeat this event each time step
    periodicity_interval = diffusion_time_step;
  }
  void step();
  void dump(const std::string indent);

  World* world;

  // this event diffuses all molecules that have this diffusion time_step
  float_t diffusion_time_step;

private:
  // molecules newly created in reactions
  std::vector<DiffuseOrUnimolReactionAction> new_diffuse_or_unimol_react_actions;

  void diffuse_molecules(Partition& p, const std::vector< molecule_id_t >& indices);

  void diffuse_single_molecule(
      Partition& p,
      const molecule_id_t vm_id,
      const float_t time_up_to_event_end
  );

  // ---------------------------------- volume molecules ----------------------------------

  void diffuse_vol_molecule(
      Partition& p,
      const molecule_id_t vm_id,
      const float_t remaining_time_step,
      bool& was_defunct,
      vec3_t& new_pos,
      subpart_index_t& new_subpart_index
  );
  bool collide_and_react_with_vol_mol(
      Partition& p,
      Collision& collision,
      vec3_t& displacement,
      float_t remaining_time_step,
      float_t r_rate_factor
  );

  int collide_and_react_with_surf_mol(
      Partition& p,
      Collision& collision,
      float_t remaining_time_step,
      float_t r_rate_factor,
      float_t current_molecule_time
  );

  // ---------------------------------- surface molecules ----------------------------------

  void diffuse_surf_molecule(
      Partition& p,
      const molecule_id_t sm_id,
      const float_t current_time,
      const float_t remaining_time_step,
      bool& was_defunct,
      vec2_t& new_loc,
      wall_index_t& new_wall_index,
      float_t& advance_time
  );

  wall_index_t ray_trace_surf(
      Partition& p,
      const Species& species,
      Molecule& sm,
      vec2_t& remaining_displacement,
      vec2_t& new_pos,
      bool& was_defunct
  );

  bool react_2D_all_neighbors(
      Partition& p,
      Molecule& sm,
      const float_t current_time,
      const float_t remaining_time_step
  );

  // ---------------------------------- reactions ----------------------------------

  int outcome_bimolecular(
      Partition& p,
      Collision& collision,
      int path,
      float_t remaining_time_step
  );

  int outcome_unimolecular(
      Partition& p,
      Molecule& vm,
      const float_t scheduled_time,
      const Reaction* unimol_rx
  );

  int outcome_products_random(
      Partition& p,
      Collision& collision,
      float_t remaining_time_step,
      int path
  );

  void create_unimol_rx_action(
      Partition& p,
      Molecule& vm,
      float_t remaining_time_step
  );

  void react_unimol_single_molecule(
      Partition& p,
      const molecule_id_t vm_id,
      const float_t scheduled_time,
      const Reaction* unimol_rx
  );
};


//FIXME: move to some utility namespace
RayTraceState ray_trace_vol(
    Partition& p,
    rng_state& rng,
    Molecule& vm, // molecule that we are diffusing, we are changing its pos  and possibly also subvolume
    const wall_index_t previous_reflected_wall, // is WALL_INDEX_INVALID when our molecule did not replect from anything this iddfusion step yet
    vec3_t& remaining_displacement, // in/out - recomputed if there was a reflection
    collision_vector_t& molecule_collisions, // possible reactions in this part of way marching, ordered by time
    vec3_t& new_position,
    subpart_index_t& new_subpartition_index
);

void sort_collisions_by_time(collision_vector_t& molecule_collisions);


} // namespace mcell

#endif // SRC4_DIFFUSE_REACT_EVENT_H_
