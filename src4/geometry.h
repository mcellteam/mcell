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

#ifndef SRC4_GEOMETRY_H_
#define SRC4_GEOMETRY_H_

#include <vector>
#include <set>

#include "defines.h"
#include "molecule.h"
#include "simulation_config.h"
#include "dyn_vertex_structs.h"
#include "wall.h"

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
      default_color(DEFAULT_COLOR),
      is_used_in_mol_rxn_counts(false)
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
  // might be set to COUNTED_VOLUME_INDEX_INTERSECTS if the direct parent of
  // this object intersects with another object
  counted_volume_index_t counted_volume_index_outside;

  // set in initialize_is_fully_transparent, true if all walls have only surface class
  // that makes them transparent to all molecules, used in overlapping wall detections
  bool is_fully_transparent;

  // default color id of this object
  rgba_t default_color;

  // if wall index is not present in this map the wall's color is the default_color
  std::map<wall_index_t, rgba_t> wall_specific_colors;

  bool represents_compartment() const {
    return vol_compartment_id != BNG::COMPARTMENT_ID_NONE;
  }

  bool is_counted_volume_or_compartment() const {
    assert(vol_compartment_id != BNG::COMPARTMENT_ID_INVALID);
    return is_used_in_mol_rxn_counts || represents_compartment();
  }

  void initialize_neighboring_walls_and_their_edges(Partition& p);

  void set_is_used_in_mol_rxn_counts(const bool value = true) {
    is_used_in_mol_rxn_counts = value;
  }

  // returns indices of all vertices in this object that are connected through an edge
  // uses caching because initial discovery is expensive
  const std::set<vertex_index_t>& get_connected_vertices(
      const Partition& p, const vertex_index_t vi);

  // returns wall index for a pair of vertex indices
  // uses caching because initial discovery is expensive
  wall_index_t get_wall_for_vertex_pair(
      const Partition& p, const vertex_index_t vi1, const vertex_index_t vi2);

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

  // cache with information on which vertices are connected through edges,
  // used in dynamic geometry collision detection
  std::map<vertex_index_t, std::set<vertex_index_t>> connected_vertices_cache;

  // cache with information on which vertex pairs belong to which walls
  std::map<UnorderedPair<vertex_index_t>, wall_index_t> vertex_pair_to_wall_cache;
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
  std::vector<InitialSurfaceReleases> initial_region_molecules;

private:
  bool volume_info_initialized;
  bool region_is_manifold;
  bool walls_per_subpart_initialized;
  bool region_waypoints_initialized;

  Vec3 bounding_box_llf;
  Vec3 bounding_box_urb;
  pos_t volume;

  WallsPerSubpartMap walls_per_subpart;

  // points known to be inside of this region, optimization
  // for checking whether a point is inside of this region
  std::set<IVec3> waypoints_in_this_region;
public:

  void initialize_volume_info_if_needed(const Partition& p);

  void compute_bounding_box(const Partition& p, Vec3& llf, Vec3& urb);

  bool is_manifold() const {
    assert(volume_info_initialized);
    return region_is_manifold;
  }

  pos_t get_volume() const {
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

typedef std::vector<GeometryObject> GeometryObjectVector;

// several utility functions related to geometry
namespace Geometry {

double compute_geometry_object_area(const Partition& p, const GeometryObject& obj);

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
