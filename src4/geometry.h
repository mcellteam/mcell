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


/* Used to transform coordinates of surface molecules diffusing between
 * adjacent walls */
class edge_t {
public:
  wall_index_t forward_index;  /* For which wall is this a forwards transform? */
  wall_index_t backward_index; /* For which wall is this a reverse transform? */

  vec2_t translate; /* Translation vector between coordinate systems */
  float_t cos_theta;         /* Cosine of angle between coordinate systems */
  float_t sin_theta;         /* Sine of angle between coordinate systems */

  //float_t length;   /* Length of the shared edge */
  //float_t length_rcp; /* Reciprocal of length of shared edge */
};


/**
 * Single instance of a wall.
 * Owned by partition, also its vertices are owned by partition.
 *
 * This is in fact a triangle, but we are keeping the naming consistent with MCell 3.
 */
class wall_t {
public:
  wall_t()
    : id(WALL_ID_INVALID), index(WALL_INDEX_INVALID), side(0), object_id(GEOMETRY_OBJECT_ID_INVALID),
      grid_index(GRID_INDEX_INVALID),
      uv_vert1_u(POS_INVALID), uv_vert2(POS_INVALID),
      area(POS_INVALID),
      normal(POS_INVALID), unit_u(POS_INVALID), unit_v(POS_INVALID), distance_to_origin(POS_INVALID)
    {
  }

  wall_id_t id; // world-unique identifier of this wall, mainly for debugging
  wall_index_t index; // index in the partition where it is contained, must be fixed if moved
  uint32_t side; // index in its parent object, not sure if really needed

  geometry_object_id_t object_id; // index of object to which this wall belongs
  //wall_class_index_t class_index; // index of this wall's class

  grid_index_t grid_index; // if is GRID_INDEX_INVALID then this wall does not have a grid, index to partition's grids

  // indices of the three triangle's vertices,
  // they are shared in a partition and a single vertex should be usually represented by just one item
  // so when a position of one vertex changes, it should affect all the triangles tht use it
  vertex_index_t vertex_indices[VERTICES_IN_TRIANGLE]; // order is important since is specifies orientation

  edge_t edges[VERTICES_IN_TRIANGLE];

  float_t uv_vert1_u;   /* Surface u-coord of 2nd corner (v=0) */
  vec2_t uv_vert2;      /* Surface coords of third corner */

  float_t area;  /* Area of this element */
  vec3_t normal; /* Normal vector for this wall */
  vec3_t unit_u; /* U basis vector for this wall */
  vec3_t unit_v; /* V basis vector for this wall */
  float_t distance_to_origin; // distance to origin (point normal form)

  // p must be the partition that contains this object
  void dump(const partition_t& p, const std::string ind) const;

  bool has_grid() const {
    return grid_index != GRID_INDEX_INVALID;
  }
};


/**
 * Surface grid.
 * Owned by partition, associated always with a single wall.
 *
 * Contains an array of tiles.
 */
class grid_t {
public:
  void initialize(const float_t area);

  uint32_t num_tiles_along_axis; // Number of slots along each axis

  uint32_t num_tiles; // Number of tiles in effector grid (triangle: grid_size^2, rectangle: 2*grid_size^2)

  float_t inv_strip_wid; /* Reciprocal of the width of one strip */
  float_t vert2_slope;   /* Slope from vertex 0 to vertex 2 */
  float_t fullslope;     /* Slope of full width of triangle */

  void set_molecule_tile(uint32_t tile_index, molecule_id_t id) {
    assert(tile_index < molecules_per_tile.size());
    assert(molecules_per_tile[tile_index] == MOLECULE_INDEX_INVALID && "Cannot overwite a molecule that is already on tile");
    molecules_per_tile[tile_index] = id;
  }

  void reset_molecule_tile(uint32_t tile_index) {
    assert(tile_index < molecules_per_tile.size());
    assert(molecules_per_tile[tile_index] != MOLECULE_INDEX_INVALID && "Cannot reset a tile that has no molecule");
    molecules_per_tile[tile_index] = MOLECULE_INDEX_INVALID;
  }

  molecule_id_t get_molecule_on_tile(uint32_t tile_index) {
    assert(tile_index < molecules_per_tile.size());
    return molecules_per_tile[tile_index];
  }

private:
  // For now, there can be just one molecule per tile,
  // value is MOLECULE_ID_INVALID when the tile is not occupied
  // originally - sm_list
  std::vector<molecule_id_t> molecules_per_tile;
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


  static vec2_t grid2uv_random(
      const wall_t& w, const grid_t& g, const uint32_t tile_index,
      rng_state& rng
  );


  static vec3_t uv2xyz(const vec2_t& a, const wall_t& w, const vec3_t& wall_vert0) {
    return vec3_t(a.u) * w.unit_u + vec3_t(a.v) * w.unit_v + wall_vert0;
  }


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
