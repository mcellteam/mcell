/******************************************************************************
 *
 * Copyright (C) 2019-2021 by
 * The Salk Institute for Biological Studies and
 * Pittsburgh Supercomputing Center, Carnegie Mellon University
 *
 * Use of this source code is governed by an MIT-style
 * license that can be found in the LICENSE file or at
 * https://opensource.org/licenses/MIT.
 *
******************************************************************************/

#ifndef SRC4_GEOMETRY_H_
#define SRC4_GEOMETRY_H_

#include <memory>
#include <vector>
#include <set>

#include "defines.h"
#include "molecule.h"
#include "simulation_config.h"
#include "dyn_vertex_structs.h"
#include "wall.h"
#include "mesh_table.hpp"

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
      encompassing_region_id(REGION_ID_INVALID),
      vol_compartment_id(BNG::COMPARTMENT_ID_NONE),
      surf_compartment_id(BNG::COMPARTMENT_ID_NONE),
      counted_volume_index_inside(COUNTED_VOLUME_INDEX_INVALID),
      counted_volume_index_outside(COUNTED_VOLUME_INDEX_INVALID),
      is_fully_transparent(false),
      has_overlapped_walls(false),
      default_color(DEFAULT_COLOR),
      first_vertex_index(VERTEX_INDEX_INVALID),
      is_used_in_mol_rxn_counts(false),
      coordinate_facts_current(false),
      current_area(0),
      current_signed_volume(0)
    {
  }

  geometry_object_id_t id; // world-unique geometry object ID
  geometry_object_index_t index; // partition-unique geometry object index

  region_id_t encompassing_region_id; // ID of Region that represents this whole object, used only in pymcell4 for now

  // ID of compartment that represents volume enclosed by this object,none by default
  BNG::compartment_id_t vol_compartment_id;
  // ID of compartment that represents the surface of this object, none by default
  BNG::compartment_id_t surf_compartment_id;

  std::string name;
  std::string parent_name;

  // all walls (triangles) that form this object,
  // wall indices are unique in each partition, i.e. they are not counted from 0 for each object
  std::vector<wall_index_t> wall_indices;

  // all regions located on this object
  uint_set<region_index_t> regions;

  // counted volume to be set when a molecule goes inside of this object
  // might be set to COUNTED_VOLUME_INDEX_INTERSECTS if this object intersects
  // with another object
  counted_volume_index_t counted_volume_index_inside;

  // counted volume to be set when a molecule goes outside of this object
  // might be set to COUNTED_VOLUME_INDEX_INTERSECTS if this object or an object
  // that contains it intersects with another object
  counted_volume_index_t counted_volume_index_outside;

  // set in initialize_is_fully_transparent, true if all walls have only surface class
  // that makes them transparent to all molecules, used in overlapping wall detections
  bool is_fully_transparent;

  // to true set in check_for_overlapped_walls iff at least one of this object's
  // wall is overlapped, such object does not allow surface releases and cannot be a
  // BNGL compartment
  bool has_overlapped_walls;

  // default color id of this object
  rgba_t default_color;

  // if wall index is not present in this map the wall's color is the default_color
  std::map<wall_index_t, rgba_t> wall_specific_colors;

  // The checks of the object's mesh as the model gives it
  // (MeshUtils::check_lists or check_parsed_walls), kept for the whole run: the
  // MeshTable, in which wall j is wall_indices[j] and vertex v is the
  // partition's vertex first_vertex_index + v, the flags (manifold,
  // consistently oriented, watertight) and the shells. The converter that
  // creates the object sets it with set_mesh_check once the walls are linked.
  std::shared_ptr<const MeshUtils::MeshCheck> mesh_check;

  // the partition's index of the object's vertex 0; the converters add an
  // object's vertices one after another
  vertex_index_t first_vertex_index;

  void set_mesh_check(std::shared_ptr<const MeshUtils::MeshCheck> check);

  const MeshUtils::MeshCheck& get_mesh_check() const {
    release_assert(mesh_check != nullptr && mesh_check->table_built);
    return *mesh_check;
  }

  // the object's vertex that is the partition's vertex vi
  MeshUtils::mesh_index_t get_object_vertex_index(const vertex_index_t vi) const;

  // the MeshTable's wall that is the object's wall wi (the partition's index):
  // table wall j is wall_indices[j]
  MeshUtils::mesh_index_t get_table_wall_index(const wall_index_t wi) const;

  // the side k' of the wall across side k of the object's wall wi, where the
  // two walls share the edge, from the MeshTable; EDGE_INDEX_INVALID when the
  // edge has no second wall
  edge_index_t get_neighbor_side(const wall_index_t wi, const edge_index_t k) const;

  // The area and the signed volume of the object's walls, in internal units;
  // the volume is 0 unless the mesh is watertight. Until a vertex of the object
  // moves they are the values of the checks; after a move they are recomputed
  // at the next call from the current vertices, with the same function
  // (GeometryAnalysis::volume_and_area) over the same walls in the same order.
  double get_area(const Partition& p) const;
  double get_signed_volume(const Partition& p) const;

  // called by Geometry::update_moved_walls when a vertex of the object moves
  void set_coordinates_changed() {
    coordinate_facts_current = false;
  }

  // the walls that use the partition's vertex vi, ascending, as indices into
  // wall_indices; the lists are built from the MeshTable on first use
  MeshUtils::IndexRange get_walls_of_vertex(const vertex_index_t vi) const;

  // the vertices that share an edge with the partition's vertex vi (its
  // 1-ring), as the object's vertex indices (add first_vertex_index), in the
  // order of the MeshTable's edges; the lists are built on first use
  MeshUtils::IndexRange get_neighbors_of_vertex(const vertex_index_t vi) const;

  bool represents_compartment() const {
    return vol_compartment_id != BNG::COMPARTMENT_ID_NONE;
  }

  bool is_counted_volume_or_compartment() const {
    assert(vol_compartment_id != BNG::COMPARTMENT_ID_INVALID);
    return is_used_in_mol_rxn_counts || represents_compartment();
  }

  // Links the walls across their shared edges from the object's MeshTable, in
  // which wall j is wall_indices[j]: the neighbour across each side and the
  // edge's change of frame. A side on a boundary edge keeps no neighbour and an
  // uninitialized edge. The walls' constants must be set.
  void link_walls(Partition& p, const MeshUtils::MeshTable& table) const;

  void set_is_used_in_mol_rxn_counts(const bool value = true) {
    is_used_in_mol_rxn_counts = value;
  }

  // the first wall, in the order of wall_indices, that has both of the
  // partition's vertices vi1 and vi2; WALL_INDEX_INVALID if there is none
  wall_index_t get_wall_for_vertex_pair(const vertex_index_t vi1, const vertex_index_t vi2) const;

  // checks all walls and their regions and if all are only transparent to all molecules,
  // sets member is_fully_transparent to true
  void initialize_is_fully_transparent(Partition& p);

  // returns nonempty string on error
  std::string validate_volumetric_mesh(const Partition& p) const;

  // p must be the partition that contains this object
  void dump(const Partition& p, const std::string ind) const;
  static void dump_array(const Partition& p, const std::vector<GeometryObject>& vec);
  void to_data_model_as_geometrical_object(
      const Partition& p, const SimulationConfig& config,
      Json::Value& object,
      std::set<rgba_t>& used_colors) const;
  void to_data_model_as_model_object(const Partition& p, Json::Value& model_object) const;

  // checks only in debug mode whether the wall index belongs to this object
  rgba_t get_wall_color(const wall_index_t wi) const;

  // sets wall_specific_colors[wi] = color even if color is default_color
  // checks only in debug mode whether the wall index belongs to this object
  void set_wall_color(const wall_index_t wi, const rgba_t color);

private:
  // true if there are MolOrRxnCountEvents that use this geometry object,
  // made private to make sure that is_counted_volume() is used when checking
  // whether crossing this objects' walls should be taken into account
  bool is_used_in_mol_rxn_counts;

  // what get_area and get_signed_volume return; recomputed when not current
  mutable bool coordinate_facts_current;
  mutable double current_area;
  mutable double current_signed_volume;
  void update_coordinate_facts(const Partition& p) const;

  // the lists per vertex, built from the MeshTable on first use
  mutable std::shared_ptr<const MeshUtils::VertexLists> walls_of_vertices;
  mutable std::shared_ptr<const MeshUtils::VertexLists> neighbors_of_vertices;
};

typedef std::map<subpart_index_t, small_vector<wall_index_t>> WallsPerSubpartMap;

// this class holds information for initial release of molecules onto regions specified by
// MDL's MODIFY_SURFACE_REGIONS/MOLECULE_DENSITY or MOLECULE_NUMBER
class InitialSurfaceReleases {
public:
  InitialSurfaceReleases(
    species_id_t species_id_,
    orientation_t orientation_,
    bool const_num_not_density_,
    uint release_num_)
    : species_id(species_id_), orientation(orientation_), const_num_not_density(const_num_not_density_),
      release_num(release_num_), release_density(FLT_INVALID) {
    assert(const_num_not_density_ && "This ctor is for const_num");
  }

  InitialSurfaceReleases(
    species_id_t species_id_,
    orientation_t orientation_,
    bool const_num_not_density_,
    double release_density_)
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
  double release_density;
};


class Region {
public:
  Region()
    : id(REGION_ID_INVALID), index(REGION_INDEX_INVALID), species_id(SPECIES_ID_INVALID),
      geometry_object_id(GEOMETRY_OBJECT_ID_INVALID), name(""),
      bounding_box_initialized(false),
      walls_per_subpart_initialized(false), region_waypoints_initialized(false)
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
  std::vector<InitialSurfaceReleases> initial_region_molecules;

private:
  // the box, the walls per subpartition and the waypoints are computed from the
  // walls' positions at their first use, and again after a vertex of the region's
  // object moves (set_coordinates_changed)
  bool bounding_box_initialized;
  bool walls_per_subpart_initialized;
  bool region_waypoints_initialized;

  Vec3 bounding_box_llf;
  Vec3 bounding_box_urb;

  WallsPerSubpartMap walls_per_subpart;

  // points known to be inside of this region, optimization
  // for checking whether a point is inside of this region
  std::set<IVec3> waypoints_in_this_region;
public:

  void initialize_bounding_box_if_needed(const Partition& p);

  void compute_bounding_box(const Partition& p, Vec3& llf, Vec3& urb);

  // called by Geometry::update_moved_walls when a vertex of the region's object
  // moves; the box, the walls per subpartition and the waypoints are computed
  // again at their next use
  void set_coordinates_changed() {
    bounding_box_initialized = false;
    walls_per_subpart_initialized = false;
    region_waypoints_initialized = false;
  }

  const Vec3& get_bounding_box_llf() const {
    assert(bounding_box_initialized);
    return bounding_box_llf;
  }

  const Vec3& get_bounding_box_urb() const {
    assert(bounding_box_initialized);
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
  bool name_has_suffix_ALL() const {
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

  // The region's borders, from the MeshTable of its object: side k of a wall in
  // the region is a border when its edge has exactly one wall in the region,
  // that is, when the edge is a boundary edge of the mesh or the wall across it
  // is not in the region. The same rule for every region, ALL included; a
  // region has no borders only when no edge of it has a wall outside it, as the
  // whole of a watertight mesh. The region's walls must be set.
  void init_border_edges(const GeometryObject& obj, const MeshUtils::MeshTable& table);

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

typedef std::vector<GeometryObject> GeometryObjectVector;

// several utility functions related to geometry
namespace Geometry {

// used when creating release event
bool compute_region_expr_bounding_box(
    World* world, const RegionExprNode* region_expr,
    Vec3& llf, Vec3& urb
);

// this is the entry point called from Partition class
void update_moved_walls(
    Partition& p,
    const std::vector<VertexMoveInfo*>& scheduled_vertex_moves,
    // we can compute all the information already from scheduled_vertex_moves,
    // but the keys of the map walls_with_their_moves are the walls that we need to update
    const WallsWithTheirMovesMap& walls_with_their_moves
);

void rgba_to_components(
    const rgba_t rgba,
    double& red, double& green, double& blue, double& alpha);

} // namespace Geometry

} // namespace MCell

#endif /* SRC4_GEOMETRY_H_ */
