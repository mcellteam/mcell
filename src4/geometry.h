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

class Partition;
class World;
class RegionExprNode;

const char* const REGION_ALL_SUFFIX_W_COMMA = ",ALL";

// counted volumes are represented as a set of all counted
// geometry objects that wholly contain the volume region
class CountedVolume {
public:
  CountedVolume()
    : index(COUNTED_VOLUME_INDEX_INVALID) {
  }

  counted_volume_index_t index;
  uint_set<geometry_object_index_t> contained_in_objects;

  // for search in a map, does not compare index
  bool operator < (const CountedVolume& other) const {
    return contained_in_objects < other.contained_in_objects;
  }

  void dump() const {
    contained_in_objects.dump();
  }
};

/**
 * A single geometrical object composed of multiple walls.
 * Vertices are accessible through the wall indices.
 * Owned by partition.
 */
class GeometryObject {
public:
  GeometryObject()
    : id(GEOMETRY_OBJECT_ID_INVALID), index(GEOMETRY_OBJECT_INDEX_INVALID),
      encompassing_region_index(REGION_ID_INVALID),
      vol_compartment_id(BNG::COMPARTMENT_ID_NONE),
      surf_compartment_id(BNG::COMPARTMENT_ID_NONE),
      counted_volume_index_inside(COUNTED_VOLUME_INDEX_INVALID),
      counted_volume_index_outside(COUNTED_VOLUME_INDEX_INVALID),
      is_used_in_mol_rxn_counts(false)
    {
  }

  geometry_object_id_t id; // world-unique geometry object ID
  geometry_object_index_t index; // partition-unique geometry object index

  region_id_t encompassing_region_index; // ID of Region that represents this whole object, used only in pymcell4 for now

  // ID of compartment that represents volume enclosed by this object,none by default
  BNG::compartment_id_t vol_compartment_id;
  // ID of compartment that represents the surface of this object, none by default
  BNG::compartment_id_t surf_compartment_id;

  std::string name;
  std::string parent_name;

  // all walls (triangles) that form this object
  std::vector<wall_index_t> wall_indices;

  // all regions located on this object
  uint_set<region_index_t> regions;

  bool represents_compartment() const {
    return vol_compartment_id != BNG::COMPARTMENT_ID_NONE;
  }

  bool is_counted_volume_or_compartment() const {
    assert(vol_compartment_id != BNG::COMPARTMENT_ID_INVALID);
    return is_used_in_mol_rxn_counts || represents_compartment();
  }

  // counted volume to be set when a molecule goes inside of this object
  // might be set to COUNTED_VOLUME_INDEX_INTERSECTS if this object intersects
  // with another object
  counted_volume_index_t counted_volume_index_inside;

  // counted volume to be set when a molecule goes outside of this object
  // might be set to COUNTED_VOLUME_INDEX_INTERSECTS if the direct parent of
  // this object intersects with another object
  counted_volume_index_t counted_volume_index_outside;

  void initialize_neighboring_walls_and_their_edges(Partition& p);

  void set_is_used_in_mol_rxn_counts(const bool value = true) {
    is_used_in_mol_rxn_counts = value;
  }

  // p must be the partition that contains this object
  void dump(const Partition& p, const std::string ind) const;
  static void dump_array(const Partition& p, const std::vector<GeometryObject>& vec);
  void to_data_model_as_geometrical_object(const Partition& p, const SimulationConfig& config, Json::Value& object) const;
  void to_data_model_as_model_object(const Partition& p, Json::Value& model_object) const;
private:
  // true if there are MolOrRxnCountEvents that use this geometry object,
  // made private to make sure that is_counted_volume() is used when checking
  // whether crossing this objects' walls should be taken into account
  bool is_used_in_mol_rxn_counts;
};

typedef std::map<subpart_index_t, small_vector<wall_index_t>> WallsPerSubpartMap;

// this class holds information for initial release of molecules onto regions specified by
// MDL's MODIFY_SURFACE_REGIONS/MOLECULE_DENSITY or MOLECULE_NUMBER
// TODO: rename to initial_surface_releases?
class InitialRegionMolecules {
public:
  InitialRegionMolecules(
    species_id_t species_id_,
    orientation_t orientation_,
    bool const_num_not_density_,
    uint release_num_)
    : species_id(species_id_), orientation(orientation_), const_num_not_density(const_num_not_density_),
      release_num(release_num_), release_density(FLT_INVALID) {
    assert(const_num_not_density_ && "This ctor is for const_num");
  }

  InitialRegionMolecules(
    species_id_t species_id_,
    orientation_t orientation_,
    bool const_num_not_density_,
    float_t release_density_)
    : species_id(species_id_), orientation(orientation_), const_num_not_density(const_num_not_density_),
      release_num(UINT_INVALID), release_density(release_density_) {
    assert(!const_num_not_density_ && "This ctor is for density");
  }

  bool is_release_by_num() const {
    return const_num_not_density;
  }

  bool is_release_by_density() const {
    return !const_num_not_density;
  }

  void to_data_model(
      const BNG::SpeciesContainer& all_species,
      Json::Value& initial_region_molecules
  ) const;

  void dump(const std::string ind) const;

  species_id_t species_id;
  orientation_t orientation;
  bool const_num_not_density;
  uint release_num;
  float_t release_density;
};


class Region {
public:
  Region()
    : id(REGION_ID_INVALID), index(REGION_INDEX_INVALID), species_id(SPECIES_ID_INVALID),
      geometry_object_id(GEOMETRY_OBJECT_ID_INVALID), name(""),
      volume_info_initialized(false), region_is_manifold(false),
      walls_per_subpart_initialized(false), region_waypoints_initialized(false),
      volume(FLT_INVALID)
      {
  }

  region_id_t id;
  region_index_t index;

  // the reactivity of a region is modeled using reactions,
  // if this region belongs to a surface class, species_id is to to represent the surface class
  // may be SPECIES_ID_INVALID is there is not surface class assigned to this region
  species_id_t species_id;

  bool has_surface_class() const {
    return species_id != SPECIES_ID_INVALID;
  }

  // to which object this region belongs
  geometry_object_id_t geometry_object_id;

  std::string name;

  // each wall contained in this map is a part of this region
  // the vector of edge indices may be empty but if not, it specifies the
  // edge od the wall that is a border of the region
  std::map<wall_index_t, std::set<edge_index_t>> walls_and_edges;

  // initial counts of molecules for this region
  std::vector<InitialRegionMolecules> initial_region_molecules;

private:
  // TODO: initialize all this when a region is created/finalized
  // and get rid of the *_initialized values

  bool volume_info_initialized;
  bool region_is_manifold;
  bool walls_per_subpart_initialized;
  bool region_waypoints_initialized;

  Vec3 bounding_box_llf;
  Vec3 bounding_box_urb;
  float_t volume;

  WallsPerSubpartMap walls_per_subpart;

  // points known to be inside of this region, optimization
  // for checking whether a point is inside of this region
  std::set<IVec3> waypoints_in_this_region;
public:

  void initialize_volume_info_if_needed(const Partition& p);

  bool is_manifold() const {
    assert(volume_info_initialized);
    return region_is_manifold;
  }

  float_t get_volume() const {
    assert(volume_info_initialized);
    return volume;
  }

  const Vec3& get_bounding_box_llf() const {
    assert(volume_info_initialized);
    return bounding_box_llf;
  }

  const Vec3& get_bounding_box_urb() const {
    assert(volume_info_initialized);
    return bounding_box_urb;
  }

  const WallsPerSubpartMap& get_walls_per_subpart() const {
    assert(walls_per_subpart_initialized);
    return walls_per_subpart;
  }

  bool is_edge(wall_index_t wall_index, edge_index_t edge_index) const {
    auto it = walls_and_edges.find(wall_index);
    if (it != walls_and_edges.end()) {
      return it->second.count(edge_index) != 0;
    }
    else {
      return false;
    }
  }

  bool is_reactive() const {
    return species_id != SPECIES_ID_INVALID;
  }

  bool is_point_inside(Partition& p, const Vec3& pos);

  // covers whole region
  // TODO: better name?
  bool name_has_all_suffix() const {
    std::string all(REGION_ALL_SUFFIX_W_COMMA);
    if (name.size() > all.size()) {
      return name.substr(name.size() - all.size()) == all;
    }
    else {
      return false;
    }
  }

  // for regions that encompass the whole geometry object ([ALL])
  void init_from_whole_geom_object(const GeometryObject& obj);

  // for regions that represent a surface region ([xxx])
  // walls must be already set
  void init_surface_region_edges(const Partition& p);

  void add_wall_to_walls_and_edges(const wall_index_t wi, const bool incl_edges = false) {
    walls_and_edges.insert(std::make_pair(wi, std::set<edge_index_t>()));
    if (incl_edges) {
      walls_and_edges[wi].insert(EDGE_INDEX_0);
      walls_and_edges[wi].insert(EDGE_INDEX_1);
      walls_and_edges[wi].insert(EDGE_INDEX_2);
    }
  }

  bool has_initial_molecules() const {
    return !initial_region_molecules.empty();
  }

  void dump(const std::string ind, const bool with_geometry = false) const;
  static void dump_array(const std::vector<Region>& vec);
  void to_data_model(const Partition& p, Json::Value& modify_surface_region) const;

private:
  void initialize_wall_subpart_mapping_if_needed(const Partition& p);

  bool initialize_region_waypoint(
      Partition& p,
      const IVec3& current_waypoint_index,
      const bool use_previous_waypoint,
      const IVec3& previous_waypoint_index,
      const bool previous_waypoint_present
  );
  void initialize_region_waypoints_if_needed(Partition& p);
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
    // FIXME: this checks should be here, but for some reason assert
    //assert(is_initialized());
    return translate;
  }

  float_t get_cos_theta() const {
    //assert(is_initialized());
    return cos_theta;
  }

  float_t get_sin_theta() const {
    //assert(is_initialized());
    return sin_theta;
  }

  void set_translate(const Vec2& value) {
    translate = value;
  }

  void set_cos_theta(const float_t value) {
    cos_theta = value;
  }

  void set_sin_theta(const float_t value) {
    sin_theta = value;
  }

private:
  // --- egde constants ---
  Vec2 translate;          /* Translation vector between coordinate systems */
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

  wall_index_t wall_index;

  uint num_tiles_along_axis; // Number of slots along each axis (originally n)
  uint num_tiles; // Number of tiles in effector grid (triangle: grid_size^2, rectangle: 2*grid_size^2) (originally n_tiles)

  float_t strip_width_rcp; /* Reciprocal of the width of one strip */ // inv_strip_wid originally
  float_t vert2_slope;   /* Slope from vertex 0 to vertex 2 */
  float_t fullslope;     /* Slope of full width of triangle */
  float_t binding_factor;
  Vec2 vert0;          /* Projection of vertex zero onto unit_u and unit_v of wall */

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

  WallCollisionRejectionData(const Vec3& normal_, const float_t& distance_to_origin_)
    : normal(normal_), distance_to_origin(distance_to_origin_) {
  }

  // for checking of consistency, the argument is usually Wall
  bool has_identical_data_as_wall(const WallCollisionRejectionData& w) const {
    return normal == w.normal && distance_to_origin == w.distance_to_origin;
  }

  Vec3 normal; /* Normal vector for this wall */
  float_t distance_to_origin; // distance to origin (point normal form)
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
      area(POS_INVALID),
      unit_u(POS_INVALID), unit_v(POS_INVALID) {

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
      area(POS_INVALID),
      unit_u(POS_INVALID), unit_v(POS_INVALID)
    {
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

  Grid grid;

  // --- wall constants ---
  bool wall_constants_precomputed;
  float_t uv_vert1_u;   /* Surface u-coord of 2nd corner (v=0) */
  Vec2 uv_vert2;      /* Surface coords of third corner */

  float_t area;  /* Area of this element */
  Vec3 unit_u; /* U basis vector for this wall */
  Vec3 unit_v; /* V basis vector for this wall */

  SubpartIndicesSet present_in_subparts; // in what subpartitions is this wall located

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
    // TODO: check for allowed value of t
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
  bool pos_is_set; // TODO: remove - covered by type
  Vec2 pos;
};


typedef std::vector<GeometryObject> GeometryObjectVector;


// several utility functions related to geometry
// TODO: move this to a separate file
namespace Geometry {


// TODO: move under Region
void compute_region_bounding_box(
    const Partition& p, const Region& r,
    Vec3& llf, Vec3& urb
);

// used when creating release event
bool compute_region_expr_bounding_box(
    World* world, const RegionExprNode* region_expr,
    Vec3& llf, Vec3& urb
);

// this is the entry point called from Partition class
void update_moved_walls(
    Partition& p,
    const std::vector<VertexMoveInfo>& scheduled_vertex_moves,
    // we can compute all the information already from scheduled_vertex_moves,
    // but the keys of the map walls_with_their_moves are the walls that we need to update
    const WallsWithTheirMovesMap& walls_with_their_moves
);

} // namespace Geometry

} // namespace MCell

#endif /* SRC4_GEOMETRY_H_ */
