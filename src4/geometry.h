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

namespace MCell {

class Partition;
class UintSet;

/**
 * A single geometrical object composed of multiple walls.
 * Vartices are accessible through the wall indices.
 * Owned by partition.
 */
class GeometryObject {
public:
  geometry_object_id_t id; // world-unique geometry object ID
  std::string name;

  // bool is_closed;

  // all walls (triangles) that form this object
  std::vector<wall_index_t> wall_indices;

  // p must be the partition that contains this object
  void dump(const Partition& p, const std::string ind) const;
};


/* Used to transform coordinates of surface molecules diffusing between
 * adjacent walls */
class Edge {
public:
  Edge()
    : forward_index(WALL_INDEX_INVALID), backward_index(WALL_INDEX_INVALID),
      edge_constants_precomputed(false), translate(0), cos_theta(0), sin_theta(0)
    {
  }

  void precompute_edge_constants(const Partition& p, int edgenum);

  wall_index_t forward_index;  /* For which wall is this a forwards transform? */
  wall_index_t backward_index; /* For which wall is this a reverse transform? */

  bool edge_constants_precomputed; // may be used only for debug

  // --- egde constants ---
  vec2_t translate;          /* Translation vector between coordinate systems */
  float_t cos_theta;         /* Cosine of angle between coordinate systems */
  float_t sin_theta;         /* Sine of angle between coordinate systems */
};


class Wall;

/**
 * Surface grid.
 * Owned by partition, associated always with a single wall.
 *
 * Contains an array of tiles.
 */
class Grid {
public:
  bool is_initialized() const {
    // Every initialized grid has at least one item in this array
    return !molecules_per_tile.empty();
  }

  void initialize(const Partition& p, const Wall& w);

  uint num_tiles_along_axis; // Number of slots along each axis (originally n)
  uint num_tiles; // Number of tiles in effector grid (triangle: grid_size^2, rectangle: 2*grid_size^2) (originally n_tiles)

  float_t strip_width_rcp; /* Reciprocal of the width of one strip */ // inv_strip_wid originally
  float_t vert2_slope;   /* Slope from vertex 0 to vertex 2 */
  float_t fullslope;     /* Slope of full width of triangle */
  float_t binding_factor;
  vec2_t vert0;          /* Projection of vertex zero onto unit_u and unit_v of wall */

  void set_molecule_tile(tile_index_t tile_index, molecule_id_t id) {
    assert(is_initialized());
    assert(num_tiles == molecules_per_tile.size());
    assert(tile_index < molecules_per_tile.size());
    assert(molecules_per_tile[tile_index] == MOLECULE_INDEX_INVALID && "Cannot overwite a molecule that is already on tile");
    molecules_per_tile[tile_index] = id;
  }

  void reset_molecule_tile(tile_index_t tile_index) {
    assert(is_initialized());
    assert(num_tiles == molecules_per_tile.size());
    assert(tile_index < molecules_per_tile.size());
    assert(molecules_per_tile[tile_index] != MOLECULE_INDEX_INVALID && "Cannot reset a tile that has no molecule");
    molecules_per_tile[tile_index] = MOLECULE_INDEX_INVALID;
  }

  molecule_id_t get_molecule_on_tile(tile_index_t tile_index) const {
    assert(is_initialized());
    assert(num_tiles == molecules_per_tile.size());
    assert(tile_index < molecules_per_tile.size());
    return molecules_per_tile[tile_index];
  }

private:
  // For now, there can be just one molecule per tile,
  // value is MOLECULE_ID_INVALID when the tile is not occupied
  // indexed by type tile_index_t
  // Every initialized grid has at least one item in this array
  std::vector<molecule_id_t> molecules_per_tile;
};


/**
 * Single instance of a wall.
 * Owned by partition, also its vertices are owned by partition.
 *
 * This is in fact a triangle, but we are keeping the naming consistent with MCell 3.
 *
 * TODO: Add additional debug checks that will make sure that the
 * state of this object is consistent. However, how to do it without
 * making the attributes private?
 */
class Wall {
public:
  Wall()
    : id(WALL_ID_INVALID), index(WALL_INDEX_INVALID), side(0), object_id(GEOMETRY_OBJECT_ID_INVALID),
      wall_constants_precomputed(false),
      uv_vert1_u(POS_INVALID), uv_vert2(POS_INVALID),
      area(POS_INVALID),
      normal(POS_INVALID), unit_u(POS_INVALID), unit_v(POS_INVALID), distance_to_origin(POS_INVALID)
    {
  }

  // the partition argument is used only to access vertices, wall is not aded to the partition
  Wall(
      const Partition& p,
      const vertex_index_t index0, const vertex_index_t index1, const vertex_index_t index2,
      const bool do_precompute_wall_constants, const bool do_precompute_edge_constants)
    : id(WALL_ID_INVALID), index(WALL_INDEX_INVALID), side(0), object_id(GEOMETRY_OBJECT_ID_INVALID),
      wall_constants_precomputed(false),
      uv_vert1_u(POS_INVALID), uv_vert2(POS_INVALID),
      area(POS_INVALID),
      normal(POS_INVALID), unit_u(POS_INVALID), unit_v(POS_INVALID), distance_to_origin(POS_INVALID)
    {
    vertex_indices[0] = index0;
    vertex_indices[1] = index1;
    vertex_indices[2] = index2;

    if (do_precompute_wall_constants) {
      precompute_wall_constants(p);
    }
    if (do_precompute_edge_constants) {
      assert(do_precompute_wall_constants);
      precompute_edge_constants(p);
    }
  }

  // needs vertex indices to be set
  void precompute_wall_constants(const Partition& p);

  // needs wall constants to be precomputed (all adjacent walls)
  void precompute_edge_constants(const Partition& p);


  wall_id_t id; // world-unique identifier of this wall, mainly for debugging
  wall_index_t index; // index in the partition where it is contained, must be fixed if moved
  uint side; // index in its parent object, not sure if really needed

  geometry_object_id_t object_id; // index of object to which this wall belongs
  //wall_class_index_t class_index; // index of this wall's class

  // indices of the three triangle's vertices,
  // they are shared in a partition and a single vertex should be usually represented by just one item
  // so when a position of one vertex changes, it should affect all the triangles that use it
  vertex_index_t vertex_indices[VERTICES_IN_TRIANGLE]; // order is important since is specifies orientation

  Edge edges[EDGES_IN_TRIANGLE]; // note: edges can be shared among walls to save memory

  // NOTE: what about walls that are neighboring over a partition edge?
  wall_index_t nb_walls[EDGES_IN_TRIANGLE]; // neighboring wall indices

  Grid grid;

  // --- wall constants ---
  bool wall_constants_precomputed;
  float_t uv_vert1_u;   /* Surface u-coord of 2nd corner (v=0) */
  vec2_t uv_vert2;      /* Surface coords of third corner */

  float_t area;  /* Area of this element */
  vec3_t normal; /* Normal vector for this wall */
  vec3_t unit_u; /* U basis vector for this wall */
  vec3_t unit_v; /* V basis vector for this wall */
  float_t distance_to_origin; // distance to origin (point normal form)

  // p must be the partition that contains this object
  void dump(const Partition& p, const std::string ind, const bool for_diff = false) const;

  bool has_initialized_grid() const {
    return grid.is_initialized();
  }

  // FIXME: since wall is contained in partition, this is not really a constant ref
  void initialize_grid(const Partition& p) {
    assert(!has_initialized_grid());
    grid.initialize(p, *this);
  }

  void update_after_vertex_change(Partition& p);
};



} /* namespace mcell */

#endif /* SRC4_GEOMETRY_H_ */
