/******************************************************************************
 *
 * Copyright (C) 2019 by
 * The Salk Institute for Biological Studies and
 * Pittsburgh Supercomputing Center, Carnegie Mellon University
 *
 * Use of this source code is governed by an MIT-style
 * license that can be found in the LICENSE file or at
 * https://opensource.org/licenses/MIT.
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

enum class CollisionType {
  INVALID,

  WALL_REDO,
  WALL_MISS,
  WALL_FRONT,
  WALL_BACK,

  VOLMOL_VOLMOL,
  SURFMOL_SURFMOL,
  VOLMOL_SURFMOL,
  UNIMOLECULAR,
  INTERMEMBRANE_SURFMOL_SURFMOL
};


class Collision;
class Partition;

// TODO: use only dense hash map/set and boost vector where possible
#ifndef INDEXER_WA
typedef boost::container::small_vector<Collision, 16> CollisionsVector;

#else
typedef std::vector<Collision> CollisionsVector;
#endif

/**
 * Information about collision of 2 volume/surface molecules or a of a wall collision,
 * used in diffuse_react and in partition.
 */
class Collision {
public:
  Collision()
#ifdef NDEBUG
    :
    type(CollisionType::INVALID), partition(nullptr),
    diffused_molecule_id(0), time(0), pos(0),
    colliding_molecule_id(0),
    rxn_class(nullptr),
    colliding_wall_index(0)
#else
    :
    type(CollisionType::INVALID), partition(nullptr),
    diffused_molecule_id(MOLECULE_ID_INVALID), time(TIME_INVALID), pos(POS_INVALID),
    colliding_molecule_id(MOLECULE_ID_INVALID),
    rxn_class(nullptr),
    colliding_wall_index(WALL_INDEX_INVALID)
#endif
    {
     }

  // vol-surf or vol-mol collision
  Collision(
      const CollisionType type_,
      Partition* partition_ptr,
      const molecule_id_t diffused_molecule_id_,
      const double time_,
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

  // surf-surf collision
  Collision(
      const CollisionType type_,
      Partition* partition_ptr,
      const molecule_id_t diffused_molecule_id_,
      const double time_, // time from event start
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
    assert(
        (type == CollisionType::SURFMOL_SURFMOL ||
         type == CollisionType::INTERMEMBRANE_SURFMOL_SURFMOL) &&
        "This constructor must be used only for surfsurf collisions");
  }

  // wall collision
  Collision(
      const CollisionType type_,
      Partition* partition_ptr,
      const molecule_id_t diffused_molecule_id_,
      const double time_,
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

  // unimolecular rxn
  Collision(
      const CollisionType type_,
      Partition* partition_ptr,
      const molecule_id_t diffused_molecule_id_,
      const double time_,
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
    assert(type == CollisionType::UNIMOLECULAR && "This constructor must be used only for unimol volmol collisions");
  }


  CollisionType type;
  Partition* partition;
  molecule_id_t diffused_molecule_id;
  double time;
  Vec3 pos;

  // valid only for is_wall_collision
  molecule_id_t colliding_molecule_id;

  // used for VOLMOL_VOLMOL or type == CollisionType::VOLMOL_SURFMOL
  BNG::RxnClass* rxn_class;

  // valid only for COLLISION_WALL*
  wall_index_t colliding_wall_index;

  bool has_pos() const {
    return
        type != CollisionType::INVALID &&
        type != CollisionType::SURFMOL_SURFMOL &&
        type != CollisionType::INTERMEMBRANE_SURFMOL_SURFMOL;
  }

  bool is_vol_mol_vol_mol_collision() const {
    return type == CollisionType::VOLMOL_VOLMOL;
  }

  bool is_unimol_reaction() const {
    return type == CollisionType::UNIMOLECULAR;
  }

  bool is_mol_mol_reaction() const {
    return
        type == CollisionType::VOLMOL_VOLMOL ||
        type == CollisionType::SURFMOL_SURFMOL ||
        type == CollisionType::VOLMOL_SURFMOL ||
        type == CollisionType::INTERMEMBRANE_SURFMOL_SURFMOL;
  }

  bool is_wall_collision() const {
    assert(type != CollisionType::WALL_REDO && "Not sure yet what to do with redo");
    return type == CollisionType::WALL_FRONT || type == CollisionType::WALL_BACK;
  }

  orientation_t get_orientation_against_wall() const {
    if (type == CollisionType::WALL_BACK) {
      return ORIENTATION_UP;
    }
    else if (type == CollisionType::WALL_FRONT) {
      return ORIENTATION_DOWN;
    }
    else {
      assert(false);
      return ORIENTATION_NOT_SET;
    }
  }

  // full dump
  void dump(Partition& p, const std::string ind) const;

  // for comparison with mcell3
  void dump(
      const Partition& p,
      const std::string extra_comment,
      const uint64_t iteration,
      double time_override = TIME_INVALID
  ) const;

  std::string to_string(const Partition& p) const;
  static void dump_array(Partition& p, const CollisionsVector& vec);
};

} /* namespace MCell */

#endif /* SRC4_COLLISION_STRUCTS_H_ */
