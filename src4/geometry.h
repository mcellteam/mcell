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
#include "molecule.h"
#include "dyn_vertex_structs.h"

namespace MCell {

class Partition;
class UintSet;


class Region {
public:

  std::string name;
  //region_index_t index; // not sure if it is needed

  // the reactivity of the region is modeled using reactions and
  // this region has its species specified
  species_id_t species_id;

  // TODO: wall indices
  // we might have a region but it is not used at all now and not referenced from anywhere
};


/**
 * A single geometrical object composed of multiple walls.
 * Vertices are accessible through the wall indices.
 * Owned by partition.
 */
class GeometryObject {
public:
  geometry_object_id_t id; // world-unique geometry object ID
  geometry_object_index_t index; // partition-unique geometry object index
  std::string name;

  // bool is_closed;

  // all walls (triangles) that form this object
  std::vector<wall_index_t> wall_indices;

  // p must be the partition that contains this object
  void dump(const Partition& p, const std::string ind) const;
};


/* Used to transform coordinates of surface molecules diffusing between
 * adjacent walls, owned by its wall */
class Edge {
public:
  Edge()
    : forward_index(WALL_INDEX_INVALID), backward_index(WALL_INDEX_INVALID),
      edge_num_used_for_init(EDGE_INDEX_INVALID), translate(0), cos_theta(0), sin_theta(0)
    {
  }

  bool is_initialized() {
    return edge_num_used_for_init != EDGE_INDEX_INVALID;
  }

  void reinit_edge_constants(const Partition& p);

  void dump() const;

  void debug_check_values_are_uptodate(const Partition& p);

  wall_index_t forward_index;  /* For which wall is this a forwards transform? */
  wall_index_t backward_index; /* For which wall is this a reverse transform? */

  // used only for debug, checks that the precomputed values are correct
  edge_index_t edge_num_used_for_init;

  // --- egde constants ---
  vec2_t translate;          /* Translation vector between coordinate systems */
  float_t cos_theta;         /* Cosine of angle between coordinate systems */
  float_t sin_theta;         /* Sine of angle between coordinate systems */
};


class Wall;

/**
 * Surface grid.
 * Owned by its wall.
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
    assert(tile_index != TILE_INDEX_INVALID);
    assert(num_tiles == molecules_per_tile.size());
    assert(tile_index < molecules_per_tile.size());
    assert(molecules_per_tile[tile_index] == MOLECULE_INDEX_INVALID && "Cannot overwite a molecule that is already on tile");

    molecules_per_tile[tile_index] = id;
    num_occupied++;
  }

  void reset_molecule_tile(tile_index_t tile_index) {
    assert(is_initialized());
    assert(tile_index != TILE_INDEX_INVALID);
    assert(num_tiles == molecules_per_tile.size());
    assert(tile_index < molecules_per_tile.size());
    assert(molecules_per_tile[tile_index] != MOLECULE_INDEX_INVALID && "Cannot reset a tile that has no molecule");

    molecules_per_tile[tile_index] = MOLECULE_INDEX_INVALID;
    num_occupied--;
  }

  molecule_id_t get_molecule_on_tile(tile_index_t tile_index) const {
    assert(is_initialized());
    assert(tile_index != TILE_INDEX_INVALID);
    assert(num_tiles == molecules_per_tile.size());
    assert(tile_index < molecules_per_tile.size());
    return molecules_per_tile[tile_index];
  }

  // populates array molecules with ids of molecules belonging to this grid
  void get_contained_molecules(
      small_vector<molecule_id_t>& molecule_ids
  ) const;

  void reset_all_tiles() {
    std::fill(molecules_per_tile.begin(), molecules_per_tile.end(), MOLECULE_INDEX_INVALID);
    num_occupied = 0;
  }

  bool is_full() const {
    assert(is_initialized());
    assert(num_occupied <= num_tiles);
    return num_occupied == num_tiles;
  }

  void dump() const;

private:
  uint num_occupied; // How many tiles are occupied

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
    : id(WALL_ID_INVALID), index(WALL_INDEX_INVALID), side(0),
      object_id(GEOMETRY_OBJECT_ID_INVALID), object_index(GEOMETRY_OBJECT_INDEX_INVALID),
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
    : id(WALL_ID_INVALID), index(WALL_INDEX_INVALID), side(0),
      object_id(GEOMETRY_OBJECT_ID_INVALID), object_index(GEOMETRY_OBJECT_INDEX_INVALID),
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
      reinit_edge_constants(p);
    }
  }

  // needs vertex indices to be set
  void precompute_wall_constants(const Partition& p);

  // needs wall constants to be precomputed (all adjacent walls)
  void reinit_edge_constants(const Partition& p);

  bool is_same_region(const Wall& w) const {
    if (w.object_id == object_id) {
      // TODO: regions
      return true;
    }
    else {
      return false;
    }
  }

  wall_id_t id; // world-unique identifier of this wall, mainly for debugging
  wall_index_t index; // index in the partition where it is contained, must be fixed if moved
  uint side; // index in its parent object, not sure if really needed

  geometry_object_id_t object_id; // id of object to which this wall belongs
  geometry_object_index_t object_index; // index of object in this partition to which this wall belongs


  // regions, one wall can belong to multiple regions, regions are owned by partition
  // TODO: not used
  std::vector<species_id_t> surface_class_species; // this is the only information that mcell keeps along this wall and it is derived from region


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

  void initialize_grid(const Partition& p) {
    assert(!has_initialized_grid());
    grid.initialize(p, *this);
  }

  void update_after_vertex_change(Partition& p);
};


// auxiliary class used to hold information on
// grid position for reactions
class GridPos {
public:
  GridPos()
    : initialized(false),
      wall_index(WALL_INDEX_INVALID), tile_index(TILE_INDEX_INVALID),
      pos_is_set(false), pos(POS_INVALID) {
  }

  static GridPos make_with_pos(const Partition& p, const Molecule& sm) {
    assert(sm.is_surf());
    assert(sm.s.pos.is_valid());
    GridPos res;

    res.initialized = true;
    res.wall_index = sm.s.wall_index;
    res.tile_index = sm.s.grid_tile_index;
    res.pos = sm.s.pos;
    res.pos_is_set = true;
    return res;
  }

  static GridPos make_without_pos(const Partition& p, const WallTileIndexPair& wall_tile_index_pair) {
    assert(wall_tile_index_pair.tile_index != TILE_INDEX_INVALID);
    GridPos res;

    res.initialized = true;
    res.wall_index = wall_tile_index_pair.wall_index;
    res.tile_index = wall_tile_index_pair.tile_index;
    res.pos_is_set = false;
    return res;
  }

  bool initialized; // was this info initialized
  wall_index_t wall_index;  /* wall where the tile is on */
  tile_index_t tile_index;  /* index on that tile */
  bool pos_is_set;
  vec2_t pos;
};

// several utility functions related to geometry
namespace Geometry {

// this is the entry point called from Partition class
void update_moved_walls(
    Partition& p,
    const std::vector<VertexMoveInfo>& scheduled_vertex_moves,
    // we can compute all the information already from scheduled_vertex_moves,
    // but the keys of the map walls_with_their_moves are the walls that we need to update
    const WallsWithTheirMovesMap& walls_with_their_moves
);

}
} /* namespace mcell */

#endif /* SRC4_GEOMETRY_H_ */
