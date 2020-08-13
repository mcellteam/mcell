/*
 * internal_structs.h
 *
 *  Created on: Aug 13, 2020
 *      Author: ahusar
 */

#ifndef LIBMCELL_API_SHARED_STRUCTS_H_
#define LIBMCELL_API_SHARED_STRUCTS_H_

#include "defines.h"

namespace MCell {

struct VertexMoveInfo {
  VertexMoveInfo(
      partition_id_t partition_id_,
      const vertex_index_t vertex_index_,
      const Vec3& displacement_)
    : partition_id(partition_id_),
      vertex_index(vertex_index_),
      displacement(displacement_) {
  }

  // id of partition of where do the move
  partition_id_t partition_id;

  // which vertex of the object to move
  vertex_index_t vertex_index;

  // and by how much
  Vec3 displacement;
};

}

#endif /* LIBMCELL_API_SHARED_STRUCTS_H_ */
