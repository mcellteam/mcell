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

#ifndef SRC4_GEOMETRY_H_
#define SRC4_GEOMETRY_H_

#include "defines.h"

namespace mcell {

class partition_t;
class subpartition_mask_t;

/**
 * A single geometrical object composed of multiple walls.
 * Vartices are accessible through the wall indices.
 * Owned by partition.
 */
class geometry_object_t {
public:
  geometry_object_id_t id; // world-unique geometry object ID
  std::string name;

  // bool is_closed;

  // all walls (triangles) that form this object
  std::vector<wall_index_t> wall_indices;

  // p must be the partition that contains this object
  void dump(const partition_t& p, const std::string ind) const;
};


/**
 * Single instance of a wall.
 * Owned by partition, also its vertices are owned by partition.
 *
 * This is in fact a triangle, but we are keeping the naming consistent with MCell 3.
 */
class wall_t {
public:
  wall_id_t id; // world-unique identifier of this wall, mainly for debugging
  uint32_t side; // index in its parent object, not sure if really needed

  geometry_object_index_t object_index; // index of object to which this wall belongs
  //wall_class_index_t class_index; // index of this wall's class

  // indices of the three triangle's vertices,
  // they are shared in a partition and a single vertex should be usually represented by just one item
  // so when a position of one vertex changes, it should affect all the triangles tht use it
  vertex_index_t vertex_indices[VERTICES_IN_TRIANGLE]; // order is important since is specifies orientation

  float_t uv_vert1_u;   /* Surface u-coord of 2nd corner (v=0) */
  vec2_t uv_vert2;      /* Surface coords of third corner */

  vec3_t normal; /* Normal vector for this wall */
  vec3_t unit_u; /* U basis vector for this wall */
  vec3_t unit_v; /* V basis vector for this wall */
  float_t distance_to_origin; // distance to origin (point normal form)

  // p must be the partition that contains this object
  void dump(const partition_t& p, const std::string ind) const;
};


/**
 * Auxiliary class that collects diverse geometry helper functions.
 * Does not have any data, therefore there is no _t suffix.
 */
class geometry {
public:
  static void wall_subparts_collision_test(
      const partition_t& p, const wall_t& w,
      subpart_indices_vector_t& colliding_subparts
  );

private:
  static void get_wall_bounding_box(
      const vec3_t w_vert[VERTICES_IN_TRIANGLE],
      vec3_t& llf, vec3_t& urb
  );
};

/*
class wall_hit_info_t {
public:
  wall_hit_info_t() : wall_index(WALL_INDEX_INVALID), time(TIME_INVALID) { }
  wall_hit_info_t(const wall_index_t wall_index_, const float_t time_, const vec3_t& pos_)
    : wall_index(wall_index_), time(time_), pos(pos_) { }

  wall_index_t wall_index;
  float_t time;
  vec3_t pos;
};*/

} /* namespace mcell */

#endif /* SRC4_GEOMETRY_H_ */
