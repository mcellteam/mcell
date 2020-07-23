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

#ifndef SRC4_COLLISION_STRUCTS_H_
#define SRC4_COLLISION_STRUCTS_H_

#include "defines.h"
#include <set>
#include "../libs/boost/container/small_vector.hpp"
#include "../libs/sparsehash/src/google/dense_hash_set"

#include "bng/bng.h"

namespace BNG {
class RxnClass;
}

namespace MCell {



enum class RayTraceState {
  // TODO: use UpperCase
  UNDEFINED,
  HIT_SUBPARTITION,
  RAY_TRACE_HIT_WALL,
  FINISHED
};


enum class CollisionType {
  INVALID,

  WALL_REDO,
  WALL_MISS,
  WALL_FRONT,
  WALL_BACK,

  VOLMOL_VOLMOL,
  SURFMOL_SURFMOL,
  VOLMOL_SURFMOL,
  UNIMOLECULAR_VOLMOL,
};


class Collision;
class Partition;

// TODO: use only dense hash map/set and boost vector where possible
#ifndef INDEXER_WA
typedef boost::container::small_vector<Collision, 16> collision_vector_t;

#else
typedef std::vector<Collision> collision_vector_t; // FIXME: shoudl be UpperCase
#endif

/**
 * Information about collision of 2 volume/surface molecules or a of a wall collision,
 * used in diffuse_react and in partition.
 */
class Collision {
public:
  Collision()
    : type(CollisionType::INVALID), partition(nullptr),
      diffused_molecule_id(MOLECULE_ID_INVALID), time(TIME_INVALID),
      colliding_molecule_id(MOLECULE_ID_INVALID),
      rxn_class(nullptr),
      colliding_wall_index(WALL_INDEX_INVALID) {
  }

  // maybe create some static constructors with better names
  Collision(
      const CollisionType type_,
      Partition* partition_ptr,
      const molecule_id_t diffused_molecule_id_,
      const float_t time_,
      const Vec3& pos_,
      const molecule_id_t colliding_molecule_id_,
      BNG::RxnClass* rxn_class_ptr
      )
    :
      type(type_),
      partition(partition_ptr),
      diffused_molecule_id(diffused_molecule_id_),
      time(time_),
      pos(pos_),
      colliding_molecule_id(colliding_molecule_id_),
      rxn_class(rxn_class_ptr),
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
      BNG::RxnClass* rxn_class_ptr
      )
    :
      type(type_),
      partition(partition_ptr),
      diffused_molecule_id(diffused_molecule_id_),
      time(time_),
      colliding_molecule_id(colliding_molecule_id_),
      rxn_class(rxn_class_ptr),
      colliding_wall_index(WALL_INDEX_INVALID) {
    assert(type == CollisionType::SURFMOL_SURFMOL && "This constructor must be used only for surfsurf collisions");
  }

  Collision(
      const CollisionType type_,
      Partition* partition_ptr,
      const molecule_id_t diffused_molecule_id_,
      const float_t time_,
      const Vec3& pos_,
      const wall_index_t colliding_wall_index_
      )
    :
      type(type_),
      partition(partition_ptr),
      diffused_molecule_id(diffused_molecule_id_),
      time(time_),
      pos(pos_),
      colliding_molecule_id(MOLECULE_ID_INVALID),
      rxn_class(nullptr),
      colliding_wall_index(colliding_wall_index_) {
    assert((type == CollisionType::WALL_BACK || type == CollisionType::WALL_FRONT) && "This constructor must be used only for wall collisions");
  }

  Collision(
      const CollisionType type_,
      Partition* partition_ptr,
      const molecule_id_t diffused_molecule_id_,
      const float_t time_,
      const Vec3& pos_,
      BNG::RxnClass* rxn_class_ptr
      )
    :
      type(type_),
      partition(partition_ptr),
      diffused_molecule_id(diffused_molecule_id_),
      time(time_),
      pos(pos_),
      colliding_molecule_id(MOLECULE_ID_INVALID),
      rxn_class(rxn_class_ptr),
      colliding_wall_index(WALL_INDEX_INVALID) {
    assert(type == CollisionType::UNIMOLECULAR_VOLMOL && "This constructor must be used only for unimol volmol collisions");
  }


  CollisionType type;
  Partition* partition;
  molecule_id_t diffused_molecule_id;
  float_t time;
  Vec3 pos;

  // valid only for is_wall_collision
  molecule_id_t colliding_molecule_id;

  // used for VOLMOL_VOLMOL or type == CollisionType::VOLMOL_SURFMOL
  // TODO: make them private? so that we can check access
  BNG::RxnClass* rxn_class;


  // valid only for COLLISION_WALL*
  wall_index_t colliding_wall_index;

  bool is_vol_mol_vol_mol_collision() const {
    return type == CollisionType::VOLMOL_VOLMOL;
  }

  bool is_unimol_reaction() const {
    return type == CollisionType::UNIMOLECULAR_VOLMOL;
  }

  bool is_mol_mol_reaction() const {
    return type == CollisionType::VOLMOL_VOLMOL || type == CollisionType::SURFMOL_SURFMOL || type == CollisionType::VOLMOL_SURFMOL;
  }

  bool is_wall_collision() const {
    assert(type != CollisionType::WALL_REDO && "Not sure yet what to do with redo");
    return type == CollisionType::WALL_FRONT || type == CollisionType::WALL_BACK;
  }

  // full dump
  void dump(Partition& p, const std::string ind) const;

  // for comparison with mcell3
  // FIXME: args are almost the same as for the variant above
  void dump(
      const Partition& p,
      const std::string extra_comment,
      const uint64_t iteration,
      float_t time_override = TIME_INVALID
  ) const;

  std::string to_string(const Partition& p) const;
  static void dump_array(Partition& p, const collision_vector_t& vec);
};

} /* namespace MCell */

#endif /* SRC4_COLLISION_STRUCTS_H_ */
