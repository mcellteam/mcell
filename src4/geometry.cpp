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

extern "C" {
#include "rng.h" // MCell 3
#include "isaac64.h"
#include "mcell_structs.h"
#include "logging.h"
}

#include <iostream>

#include "partition.h"
#include "geometry.h"

#include "geometry_utils.inc" // uses get_wall_bounding_box, maybe not include this file

using namespace std;

namespace mcell {


/***************************************************************************
distribute_wall:
  In: a wall 'w' that fully fits into patrition 'p'
  Out: colliding_subparts - indices of all the subparitions in a given partition
        where the wall is located
***************************************************************************/
// original name: distribute_wall
void geometry::wall_subparts_collision_test(
    const partition_t& p, const wall_t& w,
    subpart_indices_vector_t& colliding_subparts
) {
  vec3_t llf, urb; /* Bounding box for wall */
  float_t leeway = 1.0; /* Margin of error */

  vec3_t w_vert[VERTICES_IN_TRIANGLE];
  w_vert[0] = p.get_geometry_vertex(w.vertex_indices[0]);
  w_vert[1] = p.get_geometry_vertex(w.vertex_indices[1]);
  w_vert[2] = p.get_geometry_vertex(w.vertex_indices[2]);

  geom_util::get_wall_bounding_box(w_vert, llf, urb);

  // min
  if (llf.x < -leeway)
    leeway = -llf.x;
  if (llf.y < -leeway)
    leeway = -llf.y;
  if (llf.z < -leeway)
    leeway = -llf.z;

  // max
  if (urb.x > leeway)
    leeway = urb.x;
  if (urb.y > leeway)
    leeway = urb.y;
  if (urb.z > leeway)
    leeway = urb.z;
  leeway = EPS + leeway * EPS;
  if (p.get_world_constants().use_expanded_list) {
    leeway += p.get_world_constants().rx_radius_3d;
  }

  llf = llf - vec3_t(leeway);
  urb = urb + vec3_t(leeway);

  // let's assume for now that we are placing a cube with corners llf and urb,
  // fing what are the min and max parition indices
  ivec3_t min_subpart_indices, max_subpart_indices;
  p.get_subpart_3d_indices(llf, min_subpart_indices);
  p.get_subpart_3d_indices(urb, max_subpart_indices);

  // do we fit into just one subparition?
  if (max_subpart_indices == min_subpart_indices) {
    colliding_subparts.push_back(
        p.get_subpartition_index_from_3d_indices(min_subpart_indices)
    );
    return;
  }

#if 0
  // don't know yet what is this doing, somehow it is trying to find some subvolume
  // but don't know which one

  // this is needed for the following code (if enabled, fix loop afterwards)
  max_subpart_indices += ivec3_t(1); // max is incremented by 1 in each dimension

  cent.x = 0.33333333333 * (w_vert[0].x + w_vert[1].x + w_vert[2].x);
  cent.y = 0.33333333333 * (w_vert[0].y + w_vert[1].y + w_vert[2].y);
  cent.z = 0.33333333333 * (w_vert[0].z + w_vert[1].z + w_vert[2].z);


  for (i = x_min; i < x_max; i++) {
    if (cent.x < world->x_partitions[i])
      break;
  }
  for (j = y_min; j < y_max; j++) {
    if (cent.y < world->y_partitions[j])
      break;
  }
  for (k = z_min; k < z_max; k++) {
    if (cent.z < world->z_partitions[k])
      break;
  }

  h = (k - 1) +
      (world->nz_parts - 1) * ((j - 1) + (world->ny_parts - 1) * (i - 1));
  where_am_i = localize_wall(w, world->subvol[h].local_storage);
  if (where_am_i == NULL)
    return NULL;
#endif

  // simply insert to all subparts that cross the expanded cube for this wall
  // TODO: optimize, the wall is a 3D triangle, not a cube...
  for (int x = min_subpart_indices.x; x <= max_subpart_indices.x; x++) {
    for (int y = min_subpart_indices.y; y <= max_subpart_indices.y; y++) {
      for (int z = min_subpart_indices.z; z <= max_subpart_indices.z; z++) {
        colliding_subparts.push_back(
            p.get_subpart_index_from_3d_indices(x, y, z)
        );
      }
    }
  }
}


vec2_t geometry::grid2uv_random(
    const wall_t& w, const grid_t& g, const uint32_t tile_index,
    rng_state& rng
) {
  int root;
  int rootrem;
  int k, j, i;
  float_t over_n;
  float_t u_ran, v_ran;

  root = (int)(sqrt((float_t)tile_index));
  rootrem = tile_index - root * root;
  k = g.num_tiles_along_axis - root - 1;
  j = rootrem / 2;
  i = rootrem - 2 * j;

  over_n = 1.0 / (float_t)(g.num_tiles_along_axis);

  u_ran = rng_dbl(&rng);
  v_ran = 1.0 - sqrt(rng_dbl(&rng));

  vec2_t res;
  res.u =
      ((float_t)(j + i) + (1 - 2 * i) * (1.0 - v_ran) * u_ran) * over_n *
          w.uv_vert1_u +
      ((float_t)(k + i) + (1 - 2 * i) * v_ran) * over_n * w.uv_vert2.u;
  res.v =
      ((float_t)(k + i) + (1 - 2 * i) * v_ran) * over_n * w.uv_vert2.v;

  return res;
}


void grid_t::initialize(const float_t area) {


  num_tiles_along_axis = (int)ceil(sqrt(area));
  if (num_tiles_along_axis < 1) {
    num_tiles_along_axis = 1;
  }

  num_tiles = num_tiles_along_axis * num_tiles_along_axis;

  // sg->binding_factor = ((double)sg->n_tiles) / w->area; ??
  // TODO: init_grid_geometry(sg);

  molecules_per_tile.resize(num_tiles, MOLECULE_ID_INVALID);
}


void geometry_object_t::dump(const partition_t& p, const std::string ind) const {
  cout << ind << "geometry_object_t: id:" << id << ", name:" << name << "\n";
  for (wall_index_t i: wall_indices) {
    cout << ind << "  " << i << ": ";
    p.get_wall(i).dump(p, ind + "  ");
  }
}


void wall_t::dump(const partition_t& p, const std::string ind) const {
  cout << "id: " << id << ", side: " << side << ", object_id: " << object_id;

  for (uint32_t i = 0; i < VERTICES_IN_TRIANGLE; i++) {
    vertex_index_t vertex_index = vertex_indices[i];
    vec3_t pos = p.get_geometry_vertex(vertex_index);
    cout << ", vert_index: " << vertex_index << ":" << pos;
  }

  cout
    << ", normal " << normal
    << ", unit_u " << unit_u
    << ", unit_v " << unit_v
    << "\n";
}



} /* namespace mcell */
