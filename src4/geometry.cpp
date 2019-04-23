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

#include <iostream>

#include "partition.h"
#include "geometry.h"


using namespace std;

namespace mcell {

/***************************************************************************
wall_bounding_box:
  In: a wall
      vector to store one corner of the bounding box for that wall
      vector to store the opposite corner
  Out: No return value.  The vectors are set to define the smallest box
       that contains the wall.
***************************************************************************/
void geometry::get_wall_bounding_box(
    const vec3_t w_vert[VERTICES_IN_TRIANGLE],
    vec3_t& llf, vec3_t& urb
) {
  llf.x = urb.x = w_vert[0].x;
  llf.y = urb.y = w_vert[0].y;
  llf.z = urb.z = w_vert[0].z;

  if (w_vert[1].x < llf.x)
    llf.x = w_vert[1].x;
  else if (w_vert[1].x > urb.x)
    urb.x = w_vert[1].x;
  if (w_vert[2].x < llf.x)
    llf.x = w_vert[2].x;
  else if (w_vert[2].x > urb.x)
    urb.x = w_vert[2].x;

  if (w_vert[1].y < llf.y)
    llf.y = w_vert[1].y;
  else if (w_vert[1].y > urb.y)
    urb.y = w_vert[1].y;
  if (w_vert[2].y < llf.y)
    llf.y = w_vert[2].y;
  else if (w_vert[2].y > urb.y)
    urb.y = w_vert[2].y;

  if (w_vert[1].z < llf.z)
    llf.z = w_vert[1].z;
  else if (w_vert[1].z > urb.z)
    urb.z = w_vert[1].z;
  if (w_vert[2].z < llf.z)
    llf.z = w_vert[2].z;
  else if (w_vert[2].z > urb.z)
    urb.z = w_vert[2].z;
}

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

  get_wall_bounding_box(w_vert, llf, urb);

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
  if (max_subpart_indices - min_subpart_indices == ivec3_t(1)) {
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


void geometry_object_t::dump(const partition_t& p, const std::string ind) const {
  cout << ind << "geometry_object_t: id:" << id << ", name:" << name << "\n";
  for (wall_index_t i: wall_indices) {
    cout << ind << "  " << i << ": ";
    p.get_wall(i).dump(p, ind + "  ");
  }
}


void wall_t::dump(const partition_t& p, const std::string ind) const {
  cout << "id: " << id << ", side: " << side << ", object_index: " << object_index;

  for (uint32_t i = 0; i < VERTICES_IN_TRIANGLE; i++) {
    vertex_index_t index = vertex_indices[i];
    vec3_t pos = p.get_geometry_vertex(index);
    cout << ", vert_index: " << index << ":" << pos;
  }

  cout << ", normal " << normal << "\n";
}



} /* namespace mcell */
