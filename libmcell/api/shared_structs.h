/******************************************************************************
 *
 * Copyright (C) 2020-2021 by
 * The Salk Institute for Biological Studies
 *
 * Use of this source code is governed by an MIT-style
 * license that can be found in the LICENSE file or at
 * https://opensource.org/licenses/MIT.
 *
******************************************************************************/

#ifndef LIBMCELL_API_SHARED_STRUCTS_H_
#define LIBMCELL_API_SHARED_STRUCTS_H_

#include "defines.h"

namespace MCell {

struct VertexMoveInfo {
  // vertex index and partition defines the object index but
  // it is faster to have the object id available as well
  VertexMoveInfo(
      const partition_id_t partition_id_,
      const geometry_object_id_t geometry_object_id_,
      const vertex_index_t vertex_index_,
      const Vec3& displacement_)
    : partition_id(partition_id_),
      geometry_object_id(geometry_object_id_),
      vertex_index(vertex_index_),
      displacement(displacement_),
      vertex_walls_are_movable(true) {
  }

  // id of partition of where do the move
  partition_id_t partition_id;

  // id of the object to move
  geometry_object_id_t geometry_object_id;

  // which vertex of the object to move
  vertex_index_t vertex_index;

  // and by how much
  Vec3 displacement;

  // may be set to false in apply_vertex_moves_per_object when any of the walls to which
  // this vertex belongs is not movable
  bool vertex_walls_are_movable;
};


struct GeometryObjectWallUnorderedPair {
  GeometryObjectWallUnorderedPair(
      geometry_object_id_t geometry_object_id1_,
      wall_index_t wall_index1_,
      geometry_object_id_t geometry_object_id2_,
      wall_index_t wall_index2_) {

    // the pair with lower geom obj id is always the first one
    if (geometry_object_id1_ <= geometry_object_id2_) {
      geometry_object_id1 = geometry_object_id1_;
      wall_index1 = wall_index1_;
      geometry_object_id2 = geometry_object_id2_;
      wall_index2 = wall_index2_;
    }
    else {
      geometry_object_id1 = geometry_object_id2_;
      wall_index1 = wall_index2_;
      geometry_object_id2 = geometry_object_id1_;
      wall_index2 = wall_index1_;
    }
  }

  // some ordering for usage in std::set
  bool operator < (const GeometryObjectWallUnorderedPair& other) const {
    if (geometry_object_id1 != other.geometry_object_id1) {
      return geometry_object_id1 < other.geometry_object_id1;
    }
    else if (wall_index1 != other.wall_index1) {
      return wall_index1 < other.wall_index1;
    }
    else if (geometry_object_id2 != other.geometry_object_id2) {
      return geometry_object_id2 < other.geometry_object_id2;
    }
    else {
      return wall_index2 < other.wall_index2;
    }
  }

  geometry_object_id_t geometry_object_id1;
  wall_index_t wall_index1;

  geometry_object_id_t geometry_object_id2;
  wall_index_t wall_index2;
};

}

#endif /* LIBMCELL_API_SHARED_STRUCTS_H_ */
