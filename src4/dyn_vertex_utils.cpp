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

#include "dyn_vertex_utils.h"

#include "logging.h"
#include "partition.h"

namespace MCell {

namespace DynVertexUtils {

void move_vertex_and_walls(Partition& p, const std::vector<wall_index_t>& wall_indices, const vertex_move_info_t& move_info) {

  // remember info on original walls
  // TODO - will I need it?

  // move all molecules whose position will change by moving the wall
  // TODO - shrink test

  // move the vertex - this will also change the data that the walls are using
  vec3_t vertex_ref = p.get_geometry_vertex(move_info.vertex_index);
  vertex_ref = vertex_ref + move_info.translation_vec;
  if (! p.in_this_partition(vertex_ref) ) {
    mcell_log("Error: Crossing partitions is not supported yet.\n");
    exit(1);
  }

  // update each wall
  for (wall_index_t wall_index: wall_indices) {
    Wall& w = p.get_wall(wall_index);

    // debug check that grid does not have any surface molecules is present
    // in the method
    w.update_after_vertex_change(p);
  }

}

// this is the entry point called from Partition class
void move_vertices(Partition& p, const std::vector<vertex_move_info_t>& scheduled_vertex_moves) {

  // TODO: this can be optimized by moving multiple vertices at once
  for (const vertex_move_info_t& move_info: scheduled_vertex_moves) {

    // get all walls that this vertex uses
    const std::vector<wall_index_t>& wall_indices = p.get_walls_using_vertex(move_info.vertex_index);


  }
}

} // namespace DynVertexUtils

} // namespace MCell
