/******************************************************************************
 *
 * Copyright (C) 2006-2017 by
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

#ifndef SRC4_WALL_UTILS_INC_
#define SRC4_WALL_UTILS_INC_

/**
 * This file is directly included into diffuse_react_event.cpp.
 * The reason why this is not a standard .cpp + .h file is to gove the compiler
 * the opportunity to inline these functions into methods of diffuse&react event.
 */
#include <vector>

#include "diffuse_react_event.h"
#include "defines.h"
#include "world.h"
#include "partition.h"
#include "geometry.h"
#include "debug_config.h"

namespace mcell {

namespace wall_util {

/***********************************************************************
walls_share_full_edge:
  In: two walls
  Out: 1 if the walls share a full edge, 0 - otherwise.
       Here by "full" we mean that the shared edge has two endpoints
       that are the vertices of both walls w1 and w2.
************************************************************************/
static bool walls_share_full_edge(const partition_t& p, const wall_t& w1, const wall_t& w2) {
  uint i, k;
  uint count = 0;

  /* count number of shared vertices between two walls */
  for (i = 0; i < VERTICES_IN_TRIANGLE; i++) {
    for (k = 0; k < VERTICES_IN_TRIANGLE; k++) {

      if (!distinguishable_vec3( p.get_wall_vertex(w1, i), p.get_wall_vertex(w2, k), EPS) ) {
        count++;
      }
    }
  }

  return count == 2;
}

/*************************************************************************
find_nbr_walls_shared_one_vertex:
   In: the origin wall
       array with information about which vertices of the origin wall
          are shared with neighbor wall (they are indices in the
          global "world->walls_using_vertex" array).
   Out: linked list of the neighbor walls that have only one common
        vertex with the origin wall (not edge-to-edge walls, but
        vertex-to-vertex walls).
   Note: the "origin" wall is not included in the list
**************************************************************************/
static void find_nbr_walls_shared_one_vertex(
    const partition_t& p,
    const wall_t& origin_wall,
    vertex_index_t shared_verts[VERTICES_IN_TRIANGLE],
    wall_indices_t& neighboring_walls
)
{
  for (uint i = 0; i < VERTICES_IN_TRIANGLE; i++) {
    if (shared_verts[i] != VERTEX_INDEX_INVALID) {

      const std::vector<wall_index_t>& wall_indices = p.get_walls_using_vertex(shared_verts[i]);
      for (wall_index_t wi: wall_indices) {
        const wall_t& w = p.get_wall(wi);

        if (w.id == origin_wall.id) {
          // we do not care about current wall
          continue;
        }

        if (!walls_share_full_edge(p, origin_wall, w)) {
          neighboring_walls.push_back(w.id);
        }
      }
    }
  }
}


} // namespace wall_util

} // namespace mcell

#endif // SRC4_WALL_UTILS_INC_
