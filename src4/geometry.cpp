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

#include <iostream>

#include "bng/bng.h"

#include "rng.h" // MCell 3
#include "isaac64.h"
#include "mcell_structs_shared.h"
#include "logging.h"
#include "wall_util.h"

#include "partition.h"
#include "geometry.h"
#include "release_event.h"
#include "datamodel_defines.h"

#include "geometry_utils.h"
#include "geometry_utils.inl" // uses get_wall_bounding_box, maybe not include this file
#include "collision_utils.inl"
#include "dump_state.h"
#include "geometry_analysis.hpp"
#include "mesh_table.hpp"

using namespace std;

namespace MCell {

// the converters create an object's walls one after another
// (GeometryObject::link_walls checks it)
MeshUtils::mesh_index_t GeometryObject::get_table_wall_index(const wall_index_t wi) const {
  assert(!wall_indices.empty());
  const wall_index_t j = wi - wall_indices[0];
  release_assert(j < wall_indices.size() && wall_indices[j] == wi);
  return static_cast<MeshUtils::mesh_index_t>(j);
}


edge_index_t GeometryObject::get_neighbor_side(const wall_index_t wi, const edge_index_t k) const {
  assert(k < EDGES_IN_TRIANGLE);
  const MeshUtils::mesh_index_t across =
      get_mesh_check().table.neighbor(get_table_wall_index(wi) * EDGES_IN_TRIANGLE + k);
  if (across == MeshUtils::MESH_INDEX_NONE) {
    return EDGE_INDEX_INVALID;
  }
  return across % EDGES_IN_TRIANGLE;
}


void GeometryObject::link_walls(Partition& p, const MeshUtils::MeshTable& table) const {
  release_assert(table.n_walls() == wall_indices.size());
  for (size_t j = 0; j < wall_indices.size(); j++) {
    release_assert(wall_indices[j] == wall_indices[0] + j && "An object's walls must be created one after another");
  }

  // an edge with two walls joins the wall of its first entry, the forward wall,
  // with the wall of its second entry, the backward wall; entries are in the
  // order of the model's walls, so the first wall that uses the edge is the
  // forward one, as when walls were linked by coordinates. A boundary edge, the
  // only other kind once the checks have passed, keeps no neighbour and the
  // uninitialized edge
  for (MeshUtils::mesh_index_t e = 0; e < table.n_edges(); e++) {
    if (table.wall_count[e] != 2) {
      continue;
    }
    const MeshUtils::mesh_index_t entry0 = table.entry0[e];
    const MeshUtils::mesh_index_t entry1 = table.entry1[e];
    const edge_index_t side0 = entry0 % EDGES_IN_TRIANGLE;
    const edge_index_t side1 = entry1 % EDGES_IN_TRIANGLE;
    Wall& w0 = p.get_wall(wall_indices[entry0 / EDGES_IN_TRIANGLE]);
    Wall& w1 = p.get_wall(wall_indices[entry1 / EDGES_IN_TRIANGLE]);
    assert(w0.wall_constants_initialized && w1.wall_constants_initialized);

    w0.nb_walls[side0] = w1.index;
    w1.nb_walls[side1] = w0.index;

    Edge edge;
    edge.forward_index = w0.index;
    edge.backward_index = w1.index;
    edge.edge_num_used_for_init = side0;
    edge.reinit_edge_constants(p);
    w0.edges[side0] = edge;
    w1.edges[side1] = edge;
  }
}


wall_index_t GeometryObject::get_wall_for_vertex_pair(const vertex_index_t vi1, const vertex_index_t vi2) const {
  assert(vi1 != vi2);

  // the walls of vi1 are ascending, so the first one with vi2 is the first in
  // wall order
  const MeshUtils::MeshTable& table = get_mesh_check().table;
  const MeshUtils::mesh_index_t v2 = get_object_vertex_index(vi2);
  for (const MeshUtils::mesh_index_t j: get_walls_of_vertex(vi1)) {
    for (MeshUtils::mesh_index_t k = 0; k < VERTICES_IN_TRIANGLE; k++) {
      if (table.vertex(j * VERTICES_IN_TRIANGLE + k) == v2) {
        return wall_indices[j];
      }
    }
  }
  return WALL_INDEX_INVALID;
}


void GeometryObject::set_mesh_check(std::shared_ptr<const MeshUtils::MeshCheck> check) {
  release_assert(check != nullptr && check->table_built);
  release_assert(check->table.n_walls() == wall_indices.size());
  mesh_check = std::move(check);
  // the checks computed them from the vertices as the converter added them
  current_area = mesh_check->area;
  current_signed_volume = mesh_check->signed_volume;
  coordinate_facts_current = true;
}


MeshUtils::mesh_index_t GeometryObject::get_object_vertex_index(const vertex_index_t vi) const {
  release_assert(
      first_vertex_index != VERTEX_INDEX_INVALID && vi >= first_vertex_index &&
      vi - first_vertex_index < get_mesh_check().n_vertices);
  return static_cast<MeshUtils::mesh_index_t>(vi - first_vertex_index);
}


void GeometryObject::get_mesh_as_given(
    const Partition& p,
    vector<GeometryAnalysis::Point>& points, vector<GeometryAnalysis::Triangle>& triangles) const {

  const MeshUtils::MeshCheck& check = get_mesh_check();

  points.resize(check.n_vertices);
  for (size_t v = 0; v < check.n_vertices; v++) {
    const Vec3& pos = p.get_geometry_vertex(first_vertex_index + v);
    points[v] = {pos.x, pos.y, pos.z};
  }
  triangles.resize(wall_indices.size());
  for (size_t j = 0; j < wall_indices.size(); j++) {
    for (uint k = 0; k < VERTICES_IN_TRIANGLE; k++) {
      triangles[j][k] = static_cast<int>(check.table.vertex(static_cast<MeshUtils::mesh_index_t>(3 * j + k)));
    }
  }
}


// as MeshUtils::check_mesh computes them: every vertex of the object in order,
// and the walls in order with their vertices from the MeshTable
void GeometryObject::update_coordinate_facts(const Partition& p) const {
  vector<GeometryAnalysis::Point> points;
  vector<GeometryAnalysis::Triangle> triangles;
  get_mesh_as_given(p, points, triangles);

  const GeometryAnalysis::VolumeArea va = GeometryAnalysis::volume_and_area(points, triangles);
  current_area = va.area;
  current_signed_volume = get_mesh_check().watertight() ? va.volume : 0.0;
  coordinate_facts_current = true;
}


double GeometryObject::get_area(const Partition& p) const {
  if (!coordinate_facts_current) {
    update_coordinate_facts(p);
  }
  return current_area;
}


double GeometryObject::get_signed_volume(const Partition& p) const {
  if (!coordinate_facts_current) {
    update_coordinate_facts(p);
  }
  return current_signed_volume;
}


MeshUtils::IndexRange GeometryObject::get_walls_of_vertex(const vertex_index_t vi) const {
  const MeshUtils::mesh_index_t v = get_object_vertex_index(vi);
  if (walls_of_vertices == nullptr) {
    walls_of_vertices = make_shared<const MeshUtils::VertexLists>(
        MeshUtils::walls_of_vertices(get_mesh_check().table));
  }
  return walls_of_vertices->of(v);
}


MeshUtils::IndexRange GeometryObject::get_neighbors_of_vertex(const vertex_index_t vi) const {
  const MeshUtils::mesh_index_t v = get_object_vertex_index(vi);
  if (neighbors_of_vertices == nullptr) {
    neighbors_of_vertices = make_shared<const MeshUtils::VertexLists>(
        MeshUtils::neighbors_of_vertices(get_mesh_check().table));
  }
  return neighbors_of_vertices->of(v);
}


// checks all walls and their regions and if all are only transparent to all molecules,
// sets member is_fully_transparent to true
void GeometryObject::initialize_is_fully_transparent(Partition& p) {

  is_fully_transparent = false;

  const BNG::SpeciesRxnClassesMap* all_molecules_rxns =
      p.get_all_rxns().get_bimol_rxns_for_reactant(p.get_all_species().get_all_molecules_species_id());

  if (all_molecules_rxns == nullptr) {
    // no transparent surface classes for all molecules were defined at all
    return;
  }

  set<species_id_t> transparent_surf_classes_cache;

  // for each wall
  for (wall_index_t wi: wall_indices) {
    const Wall& w = p.get_wall(wi);

    bool is_transparent = false;

    // check all regions
    for (region_index_t ri: w.regions) {
      const Region& reg = p.get_region(ri);
      if (!reg.is_reactive()) {
        continue;
      }

      if (transparent_surf_classes_cache.count(reg.species_id) != 0) {
        is_transparent = true;
        continue;
      }

      auto reactions_it = all_molecules_rxns->find(reg.species_id);
      if (reactions_it == all_molecules_rxns->end()) {
        // no reactions for this type of region -> wall is not transparent
        return;
      }

      BNG::RxnClass* rxn_class = reactions_it->second;

      if (rxn_class->is_transparent_type()) {
        transparent_surf_classes_cache.insert(reg.species_id);
        is_transparent = true;
      }
    }

    if (!is_transparent) {
      // wall has other surface classes -> object is not fully transparent
      return;
    }
  }

  // ok, all walls are fully transparent
  is_fully_transparent = true;
}


std::string GeometryObject::validate_volumetric_mesh(const Partition& p) const {

  std::stringstream res;

  for (wall_index_t wi: wall_indices) {
    const Wall& w = p.get_wall(wi);

    // check that each wall has all neighbors
    for (uint k = 0; k < EDGES_IN_TRIANGLE; k++) {
      if (w.nb_walls[k] == WALL_INDEX_INVALID) {
        res << "Wall side " << w.side << ": neighbor wall with index " << k << " not found.\n";
      }
    }

    // check that all edges were initialized
    for (uint k = 0; k < EDGES_IN_TRIANGLE; k++) {
      if (w.edges[k].edge_num_used_for_init == EDGE_INDEX_INVALID) {
        res << "Wall side " << w.side << ": edge with index " << k << " was not initialized.\n";
      }
      if (w.edges[k].forward_index == WALL_INDEX_INVALID) {
        res << "Wall side " << w.side << ": edge with index " << k << " has no forward wall index.\n";
      }
      if (w.edges[k].backward_index == WALL_INDEX_INVALID) {
        res << "Wall side " << w.side << ": edge with index " << k << " has no backward wall index.\n";
      }
    }

  }

  return res.str();
}



void GeometryObject::dump(const Partition& p, const std::string ind) const {
  cout << ind <<
      "GeometryObject: id:" << id << ", name:" << name <<
      ", compartment_id " << vol_compartment_id <<
      ", is_used_in_mol_rxn_counts " << is_used_in_mol_rxn_counts << "\n";

  for (wall_index_t i: wall_indices) {
    cout << ind << "  " << i << ": \n";
    p.get_wall(i).dump(p, ind + "    ");
  }
}


void GeometryObject::dump_array(const Partition& p, const std::vector<GeometryObject>& vec) {
  cout << "GeometryObject array: " << (vec.empty() ? "EMPTY" : "") << "\n";

  for (size_t i = 0; i < vec.size(); i++) {
    cout << i << ":\n";
    vec[i].dump(p, "  ");
  }
}


void GeometryObject::to_data_model_as_geometrical_object(
    const Partition& p, const SimulationConfig& config,
    Json::Value& object,
    std::set<rgba_t>& used_colors) const {

  object[KEY_NAME] = DMUtils::remove_obj_name_prefix(parent_name, name);

  bool first = true; // to indicate when to use a comma

  // first note which vertices this object uses
  uint_set<vertex_index_t> used_vertex_indices;
  for (wall_index_t wall_index: wall_indices) {
    const Wall& w = p.get_wall(wall_index);

    for (vertex_index_t vertex_index: w.vertex_indices) {
      used_vertex_indices.insert(vertex_index);
    }
  }

  // then generate vertices and remember mapping
  Json::Value& vertex_list = object[KEY_VERTEX_LIST];
  map<vertex_index_t, uint> map_vertex_index_to_vertex_list_index;
  uint current_index_in_vertex_list = 0;
  for (vertex_index_t i = 0; i < p.get_geometry_vertex_count(); i++) {

    if (used_vertex_indices.count(i) == 1) {
      // define mapping vertex_index -> index in vertex array
      map_vertex_index_to_vertex_list_index[i] = current_index_in_vertex_list;
      current_index_in_vertex_list++;

      // append triple x, y, z
      Vec3 pos = p.get_geometry_vertex(i) * Vec3(p.config.length_unit);
      DMUtils::append_triplet(vertex_list, pos.x, pos.y, pos.z);
    }
  }

  // element connections - they correspond to the ordering of walls
  Json::Value& element_connections = object[KEY_ELEMENT_CONNECTIONS];
  for (wall_index_t wall_index: wall_indices) {
    const Wall& w = p.get_wall(wall_index);

    Json::Value vertex_indices;
    for (uint i = 0; i < VERTICES_IN_TRIANGLE; i++) {
      assert(map_vertex_index_to_vertex_list_index.count(w.vertex_indices[i]) == 1);
      vertex_indices.append(map_vertex_index_to_vertex_list_index[w.vertex_indices[i]]);
    }

    element_connections.append(vertex_indices);
  }

  // surface regions
  insertion_ordered_set<region_index_t> region_indices;
  map<wall_index_t, uint> map_wall_index_to_order_index_in_object;
  for (size_t i = 0; i < wall_indices.size(); i++) {

    map_wall_index_to_order_index_in_object[wall_indices[i]] = i;

    const Wall& w = p.get_wall(wall_indices[i]);
    for (region_index_t region_index: w.regions) {
      region_indices.insert_ordered(region_index);
    }
  }

  vector<Json::Value> surface_regions_vec;
  for (region_index_t region_index: region_indices.get_as_vector()) {
    const Region& reg = p.get_region(region_index);

    if (reg.name_has_suffix_ALL()) {
      continue;
    }

    Json::Value surface_region;

    string name = DMUtils::get_surface_region_name(reg.name);
    surface_region[KEY_NAME] = name;

    Json::Value include_elements;

    for (auto it: reg.walls_and_edges) {
       include_elements.append(map_wall_index_to_order_index_in_object[it.first]);
    }
    surface_region[KEY_INCLUDE_ELEMENTS] = include_elements;
    surface_regions_vec.push_back(surface_region);
  }

  // we may create the define_surface_regions only when there is at least one region
  if (!surface_regions_vec.empty()) {
    Json::Value& define_surface_regions = object[KEY_DEFINE_SURFACE_REGIONS];

    for (auto& surface_region: surface_regions_vec) {
      define_surface_regions.append(surface_region);
    }
  }

  // colors
  if (default_color != DEFAULT_COLOR || !wall_specific_colors.empty()) {
    used_colors.insert(default_color);
    std::map<rgba_t, int> colors_this_object;
    colors_this_object[default_color] = 0;

    // collect all wall_specific_colors and define their ID
    for (const auto& pair_index_color: wall_specific_colors) {
      rgba_t wall_color = pair_index_color.second;
      auto it = colors_this_object.find(wall_color);
      if (it == colors_this_object.end()) {
        // new item
        int local_index = colors_this_object.size();
        colors_this_object[wall_color] = local_index;
      }
    }

    // generate material names for this object, must be sorted by their id
    vector<rgba_t> sorted_colors;
    sorted_colors.resize(colors_this_object.size());
    for (const auto& pair_color_local_index: colors_this_object) {
      sorted_colors[pair_color_local_index.second] = pair_color_local_index.first;
    }
    Json::Value& material_names = object[KEY_MATERIAL_NAMES];
    for (const rgba_t color: sorted_colors) {
      material_names.append(DMUtils::color_to_mat_name(color));
    }

    // and finally assign element_material_indices
    Json::Value& element_material_indices = object[KEY_ELEMENT_MATERIAL_INDICES];
    for (wall_index_t wi: wall_indices) {
      rgba_t color;
      auto it = wall_specific_colors.find(wi);
      if (it != wall_specific_colors.end()) {
        color = it->second;
        used_colors.insert(color);
      }
      else {
        color = default_color;
      }
      assert(colors_this_object.count(color) == 1);
      element_material_indices.append(colors_this_object[color]);
    }
  }
  else {
    // material_names must not be empty
    object[KEY_MATERIAL_NAMES].append(Json::Value(KEY_VALUE_MEMBRANE));
  }
}


void GeometryObject::to_data_model_as_model_object(
    const Partition& p, Json::Value& model_object) const {
  model_object[KEY_DESCRIPTION] = "";
  model_object[KEY_OBJECT_SOURCE] = VALUE_BLENDER;
  model_object[KEY_DYNAMIC_DISPLAY_SOURCE] = "script";
  model_object[KEY_SCRIPT_NAME] = "";

  string obj_name = DMUtils::remove_obj_name_prefix(parent_name, name);
  model_object[KEY_NAME] = obj_name;

  // set defaults that may be overwritten
  model_object[KEY_MEMBRANE_NAME] = "";
  model_object[KEY_PARENT_OBJECT] = "";

  const BNG::BNGData& bng_data = p.get_all_species().get_bng_data();

  if (represents_compartment()) {
    model_object[KEY_IS_BNGL_COMPARTMENT] = true;

    const BNG::Compartment& comp3d = bng_data.get_compartment(vol_compartment_id);
    assert(comp3d.is_3d);
    if (comp3d.name != obj_name) {
      mcell_error(
          "For data model export, name of object %s must be the same as its compartment name %s.",
          obj_name.c_str(), comp3d.name.c_str());
    }

    if (comp3d.has_parent()) {
      const BNG::Compartment& comp_parent = bng_data.get_compartment(comp3d.parent_compartment_id);
      if (!comp_parent.is_3d) {
        // parent of this 3d object is its membrane
        model_object[KEY_MEMBRANE_NAME] = comp_parent.name;

        if (comp_parent.has_parent()) {
          const BNG::Compartment& comp3d_parent =
              bng_data.get_compartment(comp_parent.parent_compartment_id);
          assert(comp3d_parent.is_3d);
          // parent of this 3d object is its membrane
          model_object[KEY_PARENT_OBJECT] = comp3d_parent.name;
        }
      }
      else {
        // parent is this 3d object, membrane has no name
        model_object[KEY_PARENT_OBJECT] = comp_parent.name;
      }
    }
  }
  else {
    model_object[KEY_IS_BNGL_COMPARTMENT] = false;
  }

  model_object[KEY_DYNAMIC] = false;
}


// checks only in debug mode whether the wall index belongs to this object
rgba_t GeometryObject::get_wall_color(const wall_index_t wi) const {
  assert(!wall_indices.empty());
  assert(wi >= wall_indices.front() && wi <= wall_indices.back());

  auto it = wall_specific_colors.find(wi);
  if (it == wall_specific_colors.end()) {
    return default_color;
  }
  else {
    return it->second;
  }
}


// checks only in debug mode whether the wall index belongs to this object
void GeometryObject::set_wall_color(const wall_index_t wi, const rgba_t color) {
  assert(!wall_indices.empty());
  assert(wi >= wall_indices.front() && wi <= wall_indices.back());

  wall_specific_colors[wi] = color;
}


void InitialSurfaceReleases::to_data_model(
    const BNG::SpeciesContainer& all_species,
    Json::Value& initial_region_molecules
) const {
  initial_region_molecules[KEY_MOLECULE] = all_species.get(species_id).name;
  initial_region_molecules[KEY_ORIENT] = DMUtils::orientation_to_str(orientation);
  if (const_num_not_density) {
    initial_region_molecules[KEY_MOLECULE_NUMBER] = to_string(release_num);
  }
  else {
    initial_region_molecules[KEY_MOLECULE_DENSITY] = f_to_str(release_density);
  }
}


void InitialSurfaceReleases::dump(const std::string ind) const {
  cout << ind <<
      "species_id: " << species_id <<
      ", orientation: " << orientation <<
      ", const_num_not_density: " << const_num_not_density <<
      ", release_num: " << release_num <<
      ", release_density: " << release_density <<
      "\n";
}


void Region::init_from_whole_geom_object(const GeometryObject& obj) {
  name = obj.name + REGION_ALL_SUFFIX_W_COMMA;
  geometry_object_id = obj.id;

  // simply add all walls
  for (const wall_index_t& wi: obj.wall_indices) {
    add_wall_to_walls_and_edges(wi, false);
  }
}


void Region::init_border_edges(const GeometryObject& obj, const MeshUtils::MeshTable& table) {
  assert(geometry_object_id == obj.id);

  for (auto& [wi, border_sides]: walls_and_edges) {
    const MeshUtils::mesh_index_t j = obj.get_table_wall_index(wi);
    for (edge_index_t k = 0; k < EDGES_IN_TRIANGLE; k++) {
      // the entry across the side when its edge has two walls; on a boundary
      // edge there is none
      const MeshUtils::mesh_index_t across = table.neighbor(j * EDGES_IN_TRIANGLE + k);
      if (across == MeshUtils::MESH_INDEX_NONE ||
          walls_and_edges.count(obj.wall_indices[across / EDGES_IN_TRIANGLE]) == 0) {
        border_sides.insert(k);
      }
    }
  }
}


// the region's bounding box, from the walls' current positions; releases by
// concentration take the volume of the region's object from the object
// (GeometryObject::get_signed_volume)
void Region::initialize_bounding_box_if_needed(const Partition& p) {
  if (bounding_box_initialized) {
    return;
  }

  if (walls_and_edges.empty()) {
    mcell_internal_error("Region '%s' has NULL wall array!", name.c_str());
  }

  compute_bounding_box(p, bounding_box_llf, bounding_box_urb);

  bounding_box_initialized = true;
}


void Region::initialize_wall_subpart_mapping_if_needed(const Partition& p) {
  if (walls_per_subpart_initialized) {
    return;
  }

  perf() << "Preparing wall subpart mapping for region " << name << "\n";

  // first initialize region wall subpart mapping, from the subpartitions the walls
  // are in now (Partition::apply_vertex_moves updates them for moved walls, and a
  // move of the region's object drops this mapping, Region::set_coordinates_changed)
  walls_per_subpart.clear();

  for (const auto& wall_and_edges_pair: walls_and_edges) {
    const Wall& w = p.get_wall(wall_and_edges_pair.first);

    for (subpart_index_t si: w.present_in_subparts) {
      walls_per_subpart[si].push_back(w.index);
    }
  }

  perf() << "  - wall subpart mapping prepared\n";

  walls_per_subpart_initialized = true;
}

// returns true if the waypoint with current_waypoint_index is in this region
// similar code as in Partition::create_waypoint
bool Region::initialize_region_waypoint(
    Partition& p,
    const IVec3& current_waypoint_index,
    const bool use_previous_waypoint,
    const IVec3& previous_waypoint_index,
    const bool previous_waypoint_present_in_region
) {

  // waypoint is always in the center of a subpartition
  // thanks to the extra margin on each side, we might be out of the partition
  if (p.is_valid_waypoint_index(current_waypoint_index)) {
    Waypoint& waypoint = p.get_waypoint(current_waypoint_index);


    bool previous_waypoint_use_ok = false;
    uint num_crossed;
    if (use_previous_waypoint) {
      const Waypoint& previous_waypoint = p.get_waypoint(previous_waypoint_index);

      bool must_redo_test = false;

      num_crossed = CollisionUtils::get_num_crossed_region_walls(
          p, waypoint.pos, previous_waypoint.pos,
          *this, must_redo_test
      );

      if (!must_redo_test) {
        previous_waypoint_use_ok = true;
      }
    }

    if (previous_waypoint_use_ok) {
      // ok, test passed safely
      if ((num_crossed % 2 == 0 && previous_waypoint_present_in_region) ||
          (num_crossed % 2 == 1 && !previous_waypoint_present_in_region)) {

        // still in the same region
        waypoints_in_this_region.insert(current_waypoint_index);
        return true;
      }
      else {
        return false;
      }
    }

    bool must_redo_test = false;
    bool inside = false;
    do {
      inside = CollisionUtils::is_point_inside_region_no_waypoints(p, waypoint.pos, *this, must_redo_test);
      if (must_redo_test) {
        // updates values referenced by waypoint
        p.move_waypoint_because_positioned_on_wall(current_waypoint_index);
      }
    } while (must_redo_test);

    if (inside) {
      waypoints_in_this_region.insert(current_waypoint_index);
      return true;
    }
  }

  return false;
}

void Region::initialize_region_waypoints_if_needed(Partition& p) {

  if (region_waypoints_initialized) {
    return;
  }

  perf() << "Initializing waypoints for region " << name << "\n";

  assert(walls_per_subpart_initialized);

  // get initial waypoint
  waypoints_in_this_region.clear();

  IVec3 llf_waypoint_index;
  subpart_index_t subpart_index = p.get_subpart_index(bounding_box_llf);
  p.get_subpart_3d_indices_from_index(subpart_index, llf_waypoint_index);

  // then compute how many waypoints in each dimension we need to check
  Vec3 region_dims = bounding_box_urb - bounding_box_llf;
  // num_waypoints need to be incremented by 2 - for each of the side regions
  IVec3 num_waypoints = region_dims / Vec3(p.config.subpart_edge_length) + Vec3(2);

  bool previous_waypoint_present = false;
  bool use_previous_waypoint = false;
  IVec3 previous_waypoint_index;

  for (int x = 0; x < num_waypoints.x; x++) {
    for (int y = 0; y < num_waypoints.y; y++) {
      for (int z = 0; z < num_waypoints.z; z++) {

        IVec3 current_waypoint_index( llf_waypoint_index + IVec3(x, y, z) );

        previous_waypoint_present = initialize_region_waypoint(
            p, current_waypoint_index,
            use_previous_waypoint, previous_waypoint_index, previous_waypoint_present
        );

        use_previous_waypoint = true;
        previous_waypoint_index = current_waypoint_index;
      }

      // starting a new line - start from scratch
      // NOTE: maybe we do not need to restart this every line
      use_previous_waypoint = false;
      previous_waypoint_present = false;
    }
  }

  region_waypoints_initialized = true;

  perf() << "  - region waypoints initialized\n";
}



bool Region::is_point_inside(Partition& p, const Vec3& pos) {

  initialize_bounding_box_if_needed(p);
  initialize_wall_subpart_mapping_if_needed(p);
  initialize_region_waypoints_if_needed(p);

  // get a waypoint close to this position
  IVec3 waypoint_index;
  subpart_index_t subpart_index = p.get_subpart_index(pos); // TODO: get indices directly
  p.get_subpart_3d_indices_from_index(subpart_index, waypoint_index);
  const Waypoint& waypoint = p.get_waypoint(waypoint_index);

  map<geometry_object_index_t, uint> num_crossed_walls_per_object;
  bool must_redo_test = false;
  uint num_crossed = CollisionUtils::get_num_crossed_region_walls(
      p, pos, waypoint.pos,
      *this, must_redo_test
  );
  if (!must_redo_test) {
    // ok, using waypoint passed
    bool waypoint_is_inside_this_region = waypoints_in_this_region.count(waypoint_index) != 0;
    if (waypoint_is_inside_this_region) {
      return num_crossed % 2 == 0;
    }
    else {
      return num_crossed % 2 == 1;
    }
  }
  else {
    // let's try to compute the containment without waypoints as a fallback
    bool inside = CollisionUtils::is_point_inside_region_no_waypoints(p, pos, *this, must_redo_test);
    release_assert(!must_redo_test);
    return inside;
  }
}


void Region::dump(const std::string ind, const bool with_geometry) const {
  cout << ind << "Region : " <<
      "name:" << name <<
      ", species_id: " << ((species_id == SPECIES_ID_INVALID) ? string("invalid") : to_string(species_id)) <<
      "\n";

  cout << ind << "  "  << "initial_region_molecules :\n";
  for (const auto& initial_mols: initial_region_molecules) {
    initial_mols.dump(ind + "    ");
  }

  if (with_geometry) {
    for (auto& wall_it: walls_and_edges) {
      cout << ind << "  " << "wall " << wall_it.first << ", region edges: {";
      for (auto& edge: wall_it.second) {
        cout << edge << ", ";
      }
      cout << "}\n";
    }
  }
}


void Region::dump_array(const std::vector<Region>& vec) {
  cout << "Region array: " << (vec.empty() ? "EMPTY" : "") << "\n";

  for (size_t i = 0; i < vec.size(); i++) {
    cout << i << ":\n";
    vec[i].dump("  ");
  }
}


void Region::to_data_model(const Partition& p, Json::Value& modify_surface_region) const {

  if (initial_region_molecules.empty()) {
    DMUtils::add_version(modify_surface_region, VER_DM_2018_01_11_1330);
  }
  else {
    DMUtils::add_version(modify_surface_region, VER_DM_2020_07_12_1600);
  }

  modify_surface_region[KEY_DESCRIPTION] = "";
  modify_surface_region[KEY_OBJECT_NAME] =
      DMUtils::remove_obj_name_prefix(p.get_geometry_object(geometry_object_id).name);

  string region_name = DMUtils::get_region_name(name);
  if (region_name == "ALL") {
    modify_surface_region[KEY_REGION_SELECTION] = VALUE_ALL;
    modify_surface_region[KEY_REGION_NAME] = "";
  }
  else {
    modify_surface_region[KEY_REGION_SELECTION] = VALUE_SEL;
    modify_surface_region[KEY_REGION_NAME] = region_name;
  }

  modify_surface_region[KEY_NAME] = ""; // don't care

  if (species_id != SPECIES_ID_INVALID) {
    modify_surface_region[KEY_SURF_CLASS_NAME] = p.get_species(species_id).name;
  }

  if (!initial_region_molecules.empty()) {
    Json::Value& initial_region_molecules_list = modify_surface_region[KEY_INITIAL_REGION_MOLECULES_LIST];
    for (const InitialSurfaceReleases& info: initial_region_molecules) {
      Json::Value initial_region_molecules_item;
      info.to_data_model(p.get_all_species(), initial_region_molecules_item);
      initial_region_molecules_list.append(initial_region_molecules_item);
    }
  }
}


namespace Geometry {

/***************************************************************************
eval_rel_region_bbox:
  In: release expression for a 3D region release
      place to store LLF corner of the bounding box for the release
      place to store URB corner
  Out: true. Bounding box is set based on release expression (based boolean
       intersection of bounding boxes for each region). Whether the objects
       of a release of volume molecules are watertight is a usage error,
       checked by MeshUtils::check_usage.
***************************************************************************/
bool compute_region_expr_bounding_box(
    World* world, const RegionExprNode* expr,
    Vec3& llf, Vec3& urb
) {
  int count_regions_flag = 1;
  Partition& p = world->get_partition(0);

  if (expr->op == RegionExprOperator::LEAF_SURFACE_REGION) {
    Region& reg = p.get_region_by_id(expr->region_id);
    reg.compute_bounding_box(p, llf, urb);
    return true;
  }
  else if (expr->op == RegionExprOperator::LEAF_GEOMETRY_OBJECT) {
    GeometryObject& go = p.get_geometry_object(expr->geometry_object_id);
    Region& reg = p.get_region_by_id(go.encompassing_region_id);
    reg.compute_bounding_box(p, llf, urb);
    return true;
  }
  else {
    Vec3 llf_left, urb_left;
    compute_region_expr_bounding_box(world, expr->left, llf_left, urb_left);

    Vec3 llf_right, urb_right;
    compute_region_expr_bounding_box(world, expr->right, llf_right, urb_right);

    llf = llf_left;
    urb = urb_left;

    if (expr->op == RegionExprOperator::UNION) {
      if (llf.x > llf_right.x) {
        llf.x = llf_right.x;
      }
      if (llf.y > llf_right.y) {
        llf.y = llf_right.y;
      }
      if (llf.z > llf_right.z) {
        llf.z = llf_right.z;
      }
      if (urb.x < urb_right.x) {
        urb.x = urb_right.x;
      }
      if (urb.y < urb_right.y) {
        urb.y = urb_right.y;
      }
      if (urb.z < urb_right.z) {
        urb.z = urb_right.z;
      }
    }
    else if (expr->op == RegionExprOperator::DIFFERENCE) {
      // for difference/subtraction the MCell3 implementation returns
      // the llf_left/urb_left
    }
    else if (expr->op == RegionExprOperator::INTERSECT) {
      if (llf.x < llf_right.x) {
        llf.x = llf_right.x;
      }
      if (llf.y < llf_right.y) {
        llf.y = llf_right.y;
      }
      if (llf.z < llf_right.z) {
        llf.z = llf_right.z;
      }
      if (urb.x > urb_right.x) {
        urb.x = urb_right.x;
      }
      if (urb.y > urb_right.y) {
        urb.y = urb_right.y;
      }
      if (urb.z > urb_right.z) {
        urb.z = urb_right.z;
      }
    }
    else {
      assert(false);
      return false;
    }
    return true;
  }
}



// this is the entry point called from Partition class
void update_moved_walls(
    Partition& p,
    const std::vector<VertexMoveInfo*>& scheduled_vertex_moves,
    // we can compute all the information already from scheduled_vertex_moves,
    // but the keys of the map walls_with_their_moves are the walls that we need to update
    const WallsWithTheirMovesMap& walls_with_their_moves
) {

  // move all vertices
  std::set<geometry_object_id_t> moved_objects;
  for (const VertexMoveInfo* move_info: scheduled_vertex_moves) {
    // ignore moves where walls are fixed
    if (!move_info->vertex_walls_are_movable) {
      continue;
    }

    Vec3& vertex_ref = p.get_geometry_vertex(move_info->vertex_index);
    vertex_ref = vertex_ref + move_info->displacement;
    if (! p.in_this_partition(vertex_ref) ) {
      mcell_log("Error: Crossing partitions is not supported yet.\n");
      exit(1);
    }

    // the object's area and volume are recomputed when next asked for
    p.get_geometry_object_by_id(move_info->geometry_object_id).set_coordinates_changed();
    moved_objects.insert(move_info->geometry_object_id);
  }

  if (!moved_objects.empty()) {
    // releases compare this count with the one their box and wall areas were
    // computed for (ReleaseEvent::update_region_data_after_vertex_moves)
    p.increment_vertex_moves_count();

    // the regions of the moved objects compute their box, walls per subpartition
    // and waypoints again at their next use
    for (region_index_t ri = 0; ri < p.get_regions().size(); ri++) {
      Region& reg = p.get_region(ri);
      if (moved_objects.count(reg.geometry_object_id) != 0) {
        reg.set_coordinates_changed();
      }
    }
  }

  // update walls
  for (auto it: walls_with_their_moves) {
    wall_index_t wall_index = it.first;
    Wall& w = p.get_wall(wall_index);

    // first we need to update all wall constants
    w.initialize_wall_constants(p);

    // then update wall_collision_rejection_data that hold a copy
    p.update_wall_collision_rejection_data(w);

    // reinitialize grid
    if (w.grid.is_initialized()) {
      w.grid.initialize(p, w);
    }
  }


  // edges need to be fixed after all wall have been moved
  // otherwise the edge initialization would be using
  // inconsistent data
  // we need to update also edges of neighboring walls
  uint_set<wall_index_t> walls_to_be_updated;
  // NOTE: we might consider sharing edges in the same way as in MCell3
  for (auto it: walls_with_their_moves) {
    wall_index_t wall_index = it.first;
    Wall& w = p.get_wall(wall_index);

    walls_to_be_updated.insert(wall_index);
    for (uint n = 0; n < EDGES_IN_TRIANGLE; n++) {
      walls_to_be_updated.insert(w.nb_walls[n]);
    }
  }
  for (wall_index_t wall_index: walls_to_be_updated) {
    if (wall_index == WALL_INDEX_INVALID) {
      continue;
    }
    Wall& w = p.get_wall(wall_index);
    w.initialize_edge_constants(p);
  }
}


void rgba_to_components(
    const rgba_t rgba,
    double& red, double& green, double& blue, double& alpha) {

  const double MAX = 255.0;
  red = (((uint)rgba >> 24) & 0xFF) / MAX;
  green = (((uint)rgba >> 16) & 0xFF) / MAX;
  blue = (((uint)rgba >> 8) & 0xFF) / MAX;
  alpha = ((uint)rgba & 0xFF) / MAX;
}

} /* namespace Geometry */

} /* namespace MCell */
