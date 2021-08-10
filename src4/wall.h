/******************************************************************************
 *
 * Copyright (C) 2019 by
 * The Salk Institute for Biological Studies and
 * Pittsburgh Supercomputing Center, Carnegie Mellon University
 *
 * Use of this source code is governed by an MIT-style
 * license that can be found in the LICENSE file or at
 * https://opensource.org/licenses/MIT.
 *
******************************************************************************/

#ifndef SRC4_WALL_H_
#define SRC4_WALL_H_

#include <vector>
#include <set>

#include "defines.h"
#include "molecule.h"
#include "simulation_config.h"
#include "dyn_vertex_structs.h"

namespace Json {
class Value;
}

namespace MCell {

/* Used to transform coordinates of surface molecules diffusing between
 * adjacent walls, owned by its wall */
class Edge {
public:
  Edge()
    : forward_index(WALL_INDEX_INVALID), backward_index(WALL_INDEX_INVALID),
      edge_num_used_for_init(EDGE_INDEX_INVALID), translate(0), cos_theta(0), sin_theta(0)
    {
  }

  bool is_initialized() const {
    return edge_num_used_for_init != EDGE_INDEX_INVALID;
  }

  bool is_shared_edge() const {
    return forward_index != WALL_INDEX_INVALID && backward_index != WALL_INDEX_INVALID;
  }

  // must not be called on non-shared edge
  void reinit_edge_constants(const Partition& p);

  void dump(const std::string ind = "") const;

  void debug_check_values_are_uptodate(const Partition& p);

  wall_index_t forward_index;  /* For which wall is this a forwards transform? */
  wall_index_t backward_index; /* For which wall is this a reverse transform? */

  // used only for debug, checks that the precomputed values are correct
  edge_index_t edge_num_used_for_init;

  const Vec2& get_translate() const {
    return translate;
  }

  pos_t get_cos_theta() const {
    return cos_theta;
  }

  pos_t get_sin_theta() const {
    return sin_theta;
  }

  void set_translate(const Vec2& value) {
    translate = value;
  }

  void set_cos_theta(const pos_t value) {
    cos_theta = value;
  }

  void set_sin_theta(const pos_t value) {
    sin_theta = value;
  }

private:
  // --- egde constants ---
  Vec2 translate;          /* Translation vector between coordinate systems */
  pos_t cos_theta;         /* Cosine of angle between coordinate systems */
  pos_t sin_theta;         /* Sine of angle between coordinate systems */
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

  wall_index_t wall_index;

  uint num_tiles_along_axis; // Number of slots along each axis (originally n)
  uint num_tiles; // Number of tiles in effector grid (triangle: grid_size^2, rectangle: 2*grid_size^2) (originally n_tiles)

  pos_t strip_width_rcp; /* Reciprocal of the width of one strip */ // inv_strip_wid originally
  pos_t vert2_slope;   /* Slope from vertex 0 to vertex 2 */
  pos_t fullslope;     /* Slope of full width of triangle */

  Vec2 vert0;          /* Projection of vertex zero onto unit_u and unit_v of wall */

  pos_t binding_factor;

  void set_molecule_tile(tile_index_t tile_index, molecule_id_t id) {
    assert(is_initialized());
    assert(tile_index != TILE_INDEX_INVALID);
    assert(num_tiles == molecules_per_tile.size());
    assert(tile_index < molecules_per_tile.size());
    assert(molecules_per_tile[tile_index] == MOLECULE_ID_INVALID && "Cannot overwite a molecule that is already on tile");

    molecules_per_tile[tile_index] = id;
    num_occupied++;
  }

  void reset_molecule_tile(tile_index_t tile_index) {
    assert(is_initialized());
    assert(tile_index != TILE_INDEX_INVALID);
    assert(num_tiles == molecules_per_tile.size());
    assert(tile_index < molecules_per_tile.size());
    assert(molecules_per_tile[tile_index] != MOLECULE_ID_INVALID && "Cannot reset a tile that has no molecule");

    molecules_per_tile[tile_index] = MOLECULE_ID_INVALID;
    num_occupied--;
  }

  // returns MOLECULE_ID_INVALID when tile is not occupied
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

  const MoleculeIdsVector& get_molecules_per_tile() const {
    return molecules_per_tile;
  }

  void reset_all_tiles() {
    std::fill(molecules_per_tile.begin(), molecules_per_tile.end(), MOLECULE_ID_INVALID);
    num_occupied = 0;
  }

  bool is_full() const {
    assert(is_initialized());
    assert(num_occupied <= num_tiles);
    return num_occupied == num_tiles;
  }

  uint get_num_occupied_tiles() const {
    return num_occupied;
  }

  uint get_num_free_tiles() const {
    assert(num_tiles >= num_occupied);
    return num_tiles - num_occupied;
  }

  void dump() const;

private:
  uint num_occupied; // How many tiles are occupied

  // For now, there can be just one molecule per tile,
  // value is MOLECULE_ID_INVALID when the tile is not occupied
  // indexed by type tile_index_t
  // Every initialized grid has at least one item in this array
  MoleculeIdsVector molecules_per_tile;
};


/**
 * Data that are used for fast rejection of wall collisions, 
 * stored as a copy in partition to optimize cache performance.
 */
class WallCollisionRejectionData {
public:
  WallCollisionRejectionData()
    : normal(POS_INVALID), distance_to_origin(POS_INVALID) {
  }

  WallCollisionRejectionData(const Vec3& normal_, const pos_t& distance_to_origin_)
    : normal(normal_), distance_to_origin(distance_to_origin_) {
  }

  // for checking of consistency, the argument is usually Wall
  bool has_identical_data_as_wall(const WallCollisionRejectionData& w) const {
    return normal == w.normal && distance_to_origin == w.distance_to_origin;
  }

  Vec3 normal; /* Normal vector for this wall */
  pos_t distance_to_origin; // distance to origin (point normal form)
};

/**
 * Single instance of a wall.
 * Owned by partition, also its vertices are owned by partition.
 *
 * This is in fact a triangle, but we are keeping the naming consistent with MCell 3.
 *
 * TODO_LATER: Add additional debug checks that will make sure that the
 * state of this object is consistent. However, how to do it without
 * making the attributes private?
 */
class Wall : public WallCollisionRejectionData {
public:
  Wall()
    : id(WALL_ID_INVALID), index(WALL_INDEX_INVALID), side(0),
      object_id(GEOMETRY_OBJECT_ID_INVALID), object_index(GEOMETRY_OBJECT_INDEX_INVALID),
      is_movable(true),
      wall_constants_precomputed(false),
      uv_vert1_u(POS_INVALID), uv_vert2(POS_INVALID),
      unit_u(POS_INVALID), unit_v(POS_INVALID),
      area(POS_INVALID) {

    nb_walls[0] = WALL_INDEX_INVALID;
    nb_walls[1] = WALL_INDEX_INVALID;
    nb_walls[2] = WALL_INDEX_INVALID;
  }

  // the partition argument is used only to access vertices, wall is not aded to the partition
  Wall(
      const Partition& p,
      const vertex_index_t index0, const vertex_index_t index1, const vertex_index_t index2,
      const bool do_precompute_wall_constants, const bool do_precompute_edge_constants)
    : id(WALL_ID_INVALID), index(WALL_INDEX_INVALID), side(0),
      object_id(GEOMETRY_OBJECT_ID_INVALID), object_index(GEOMETRY_OBJECT_INDEX_INVALID),
      is_movable(true),
      wall_constants_precomputed(false),
      uv_vert1_u(POS_INVALID), uv_vert2(POS_INVALID),
      unit_u(POS_INVALID), unit_v(POS_INVALID),
      area(POS_INVALID) {
    vertex_indices[0] = index0;
    vertex_indices[1] = index1;
    vertex_indices[2] = index2;

    nb_walls[0] = WALL_INDEX_INVALID;
    nb_walls[1] = WALL_INDEX_INVALID;
    nb_walls[2] = WALL_INDEX_INVALID;

    if (do_precompute_wall_constants) {
      precompute_wall_constants(p);
    }
    if (do_precompute_edge_constants) {
      assert(do_precompute_wall_constants);
      reinit_edge_constants(p);
    }
  }

  bool exists_in_partition() const {
    return id != WALL_ID_NOT_IN_PARTITION;
  }

  // needs vertex indices to be set
  void precompute_wall_constants(const Partition& p);

  // needs wall constants to be precomputed (all adjacent walls)
  void reinit_edge_constants(const Partition& p);

  bool is_same_region(const Wall& w) const {
    if (w.object_id == object_id && w.regions == regions) {
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

  bool is_movable;

  // regions, one wall can belong to multiple regions, regions are owned by partition
  // region may represent a surface class
  // all regions listed here must be a part of 'regions' in the parwent geometry object
  uint_set<region_index_t> regions;

  // indices of the three triangle's vertices,
  // they are shared in a partition and a single vertex should be usually represented by just one item
  // so when a position of one vertex changes, it should affect all the triangles that use it
  vertex_index_t vertex_indices[VERTICES_IN_TRIANGLE]; // order is important since is specifies orientation

  // note: edges can be shared among walls to save memory
  // also they may be stored using std::vector to make the Wall object smaller
  Edge edges[EDGES_IN_TRIANGLE];

  // NOTE: what about walls that are neighboring over a partition edge?
  wall_index_t nb_walls[EDGES_IN_TRIANGLE]; // neighboring wall indices

  SubpartIndicesSet present_in_subparts; // in what subpartitions is this wall located

  Grid grid;

  // --- wall constants ---
  bool wall_constants_precomputed;
  pos_t uv_vert1_u;   /* Surface u-coord of 2nd corner (v=0) */
  Vec2 uv_vert2;      /* Surface coords of third corner */

  Vec3 unit_u; /* U basis vector for this wall */
  Vec3 unit_v; /* V basis vector for this wall */

  pos_t area;  /* Area of this element (triangle) */

  // p must be the partition that contains this object
  void dump(const Partition& p, const std::string ind, const bool for_diff = false) const;

  bool has_initialized_grid() const {
    return grid.is_initialized();
  }

  void initialize_grid(const Partition& p) {
    assert(!has_initialized_grid());
    grid.initialize(p, *this);
  }
};


// variant of a wall but this version has its vertices stored in its object,
// not in the partition,
// a Wall reference based on WallWithVertices has its id == WALL_ID_NOT_IN_PARTITION because it was
// not inserted into a partition
class WallWithVertices: public Wall {
public:
  // edge constants initialization is not supported
  WallWithVertices(const Partition& p,
      const Vec3& v0, const Vec3& v1, const Vec3& v2,
      const bool do_precompute_wall_constants) {

    id = WALL_ID_NOT_IN_PARTITION;

    vertices[0] = v0;
    vertices[1] = v1;
    vertices[2] = v2;

    if (do_precompute_wall_constants) {
      precompute_wall_constants(p);
    }
  }

  Vec3 vertices[VERTICES_IN_TRIANGLE];
};


enum class GridPosType {
  NOT_INITIALIZED,
  NOT_ASSIGNED, // already holds information on position but not used by product
  REACA_UV,
  REACB_UV,
  POS_UV,
  //USE_REACB_UV,
  RANDOM
};

class GridPos;
typedef std::vector<GridPos> GridPosVector;

// auxiliary class used to hold information on
// grid position for reactions in DiffuseReactEvent::find_surf_product_positions
class GridPos {
public:
  GridPos()
    : type(GridPosType::NOT_INITIALIZED),
      wall_index(WALL_INDEX_INVALID), tile_index(TILE_INDEX_INVALID),
      pos_is_set(false), pos(POS_INVALID) {
  }

  static GridPos make_with_pos(const Vec2& pos, const WallTileIndexPair& wall_tile_index_pair) {
    GridPos res;
    res.type = GridPosType::NOT_ASSIGNED;
    res.wall_index = wall_tile_index_pair.wall_index;
    res.tile_index = wall_tile_index_pair.tile_index;
    res.pos = pos;
    res.pos_is_set = true;
    return res;
  }

  static GridPos make_with_mol_pos(const Molecule& sm) {
    assert(sm.is_surf());
    assert(sm.s.pos.is_valid());
    GridPos res;

    res.type = GridPosType::NOT_ASSIGNED;
    res.wall_index = sm.s.wall_index;
    res.tile_index = sm.s.grid_tile_index;
    res.pos = sm.s.pos;
    res.pos_is_set = true;
    return res;
  }

  static GridPos make_without_pos(const WallTileIndexPair& wall_tile_index_pair) {
    assert(wall_tile_index_pair.tile_index != TILE_INDEX_INVALID);
    GridPos res;

    res.type = GridPosType::NOT_ASSIGNED;
    res.wall_index = wall_tile_index_pair.wall_index;
    res.tile_index = wall_tile_index_pair.tile_index;
    res.pos_is_set = false;
    return res;
  }

  // returns null if not found, asserts in debug mode if there are multiple options
  static const GridPos* get_second_surf_product_pos(const GridPosVector& vec, const uint first_index) {
    const GridPos* res = nullptr;
    for (uint i = 0; i < vec.size(); i++) {
      if (i != first_index && vec[i].is_assigned()) {
        assert(res == nullptr);
        res = &vec[i];
      }
    }
    return res;
  }

  bool is_initialized() const {
    return type != GridPosType::NOT_INITIALIZED;
  }

  bool is_assigned() const {
    return is_initialized() && type != GridPosType::NOT_ASSIGNED;
  }

  void set_reac_type(const GridPosType t) {
    assert(is_initialized() && !is_assigned());
    assert(t >= GridPosType::REACA_UV && t <= GridPosType::RANDOM);
    type = t;
  }

  bool has_same_wall_and_grid(const Molecule& m) const {
    assert(m.is_surf());
    return m.s.wall_index == wall_index && m.s.grid_tile_index == tile_index;
  }

  bool has_same_wall_and_grid(const GridPos& gp) const {
    return gp.wall_index == wall_index && gp.tile_index == tile_index;
  }

  GridPosType type; // was this info initialized
  wall_index_t wall_index;  /* wall where the tile is on */
  tile_index_t tile_index;  /* index on that tile */
  bool pos_is_set;
  Vec2 pos;
};

} // namespace MCell

#endif /* SRC4_WALL_H_ */
