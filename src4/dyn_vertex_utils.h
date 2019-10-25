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

#ifndef SRC4_DYN_VERTEX_UTILS_H_
#define SRC4_DYN_VERTEX_UTILS_H_

#include "defines.h"
#include <math.h>

namespace MCell {


class Partition;

struct VertexMoveInfo {
  VertexMoveInfo(const vertex_index_t vertex_index_, const vec3_t& translation_vec_)
    : vertex_index(vertex_index_), translation_vec(translation_vec_) {
  }
  // which index to move
  vertex_index_t vertex_index;
  // and by how much
  vec3_t translation_vec;
};

typedef std::vector<VertexMoveInfo> VertexMoveInfoVector;

namespace DynVertexUtils {



void move_vertices(Partition& p, const std::vector<VertexMoveInfo>& scheduled_vertex_moves);

} // namespace DynVertexUtils

} // namespace MCell

#endif /* SRC4_DYN_VERTEX_UTILS_H_ */
