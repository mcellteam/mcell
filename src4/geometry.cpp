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

#include "bng/bng.h"

#include "rng.h" // MCell 3
#include "isaac64.h"
#include "mcell_structs.h"
#include "logging.h"

#include "partition.h"
#include "geometry.h"
#include "datamodel_defines.h"

#include "geometry_utils.inc" // uses get_wall_bounding_box, maybe not include this file

using namespace std;

namespace MCell {


void GeometryObject::dump(const Partition& p, const std::string ind) const {
  cout << ind << "GeometryObject: id:" << id << ", name:" << name << "\n";
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

void GeometryObject::to_data_model(
    const Partition& p, const SimulationConfig& config, Json::Value& object) const {

  vector<uint> vertex_order;

  object[KEY_NAME] = DMUtil::remove_obj_name_prefix(parent_name, name);
  object[KEY_MATERIAL_NAMES].append(Json::Value(KEY_VALUE_MEMBRANE));


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
      DMUtil::json_append_triplet(vertex_list, pos.x, pos.y, pos.z);
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

    if (reg.name_has_all_suffix()) {
      continue;
    }

    Json::Value surface_region;

    string name = DMUtil::get_surface_region_name(reg.name);
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

}

void Region::dump(const std::string ind) const {
  cout << ind << "Region : " <<
      "name:" << name <<
      ", species_id: " << ((species_id == SPECIES_ID_INVALID) ? string("invalid") : to_string(species_id)) <<
      "\n";

  for (auto& wall_it: walls_and_edges) {
    cout << ind << "  " << "wall " << wall_it.first << ", region edges: {";
    for (auto& edge: wall_it.second) {
      cout << edge << ", ";
    }
    cout << "}\n";
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
  DMUtil::json_add_version(modify_surface_region, JSON_DM_VERSION_1330);
  modify_surface_region[KEY_DESCRIPTION] = "";
  modify_surface_region[KEY_OBJECT_NAME] =
      DMUtil::remove_obj_name_prefix(p.get_geometry_object(geometry_object_id).name);

  string region_name = DMUtil::get_region_name(name);
  if (region_name == "ALL") {
    modify_surface_region[KEY_REGION_SELECTION] = VALUE_ALL;
    modify_surface_region[KEY_REGION_NAME] = "";
  }
  else {
    modify_surface_region[KEY_REGION_SELECTION] = VALUE_SEL;
    modify_surface_region[KEY_REGION_NAME] = region_name;
  }

  modify_surface_region[KEY_NAME] = ""; // don't care
  CONVERSION_CHECK(species_id != SPECIES_ID_INVALID, "Can convert only reactive regions");
  modify_surface_region[KEY_SURF_CLASS_NAME] = p.get_all_species().get(species_id).name;
}



// may be also used for reinitialization
void Grid::initialize(const Partition& p, const Wall& w) {

  if (is_initialized()) {
    // keep the same number of items
    #ifndef NDEBUG
      small_vector<molecule_id_t> molecule_ids;
      get_contained_molecules(molecule_ids);
      assert(molecule_ids.size() == num_occupied);
    #endif
  }
  else {
    num_occupied = 0;
  }

  num_tiles_along_axis = (int)ceil_f(sqrt_f(w.area));
  if (num_tiles_along_axis < 1) {
    num_tiles_along_axis = 1;
  }

  num_tiles = num_tiles_along_axis * num_tiles_along_axis;

  molecules_per_tile.resize(num_tiles, MOLECULE_ID_INVALID);

  strip_width_rcp = 1.0 / (w.uv_vert2.v / ((float_t)num_tiles_along_axis));
  vert2_slope = w.uv_vert2.u / w.uv_vert2.v;
  fullslope = w.uv_vert1_u / w.uv_vert2.v;

  binding_factor = ((float_t)num_tiles) / w.area;

  const Vec3& vert0_tmp = p.get_wall_vertex(w, 0);

  vert0.u = dot(vert0_tmp, w.unit_u);
  vert0.v = dot(vert0_tmp, w.unit_v);
}


// populates array molecules with ids of molecules belonging to this grid
void Grid::get_contained_molecules(
    small_vector<molecule_id_t>& molecule_ids
) const {
  // might be optimized by retaining a vector/set that gets invalidated if anything changes
  molecule_ids.clear();
  for (molecule_id_t id: molecules_per_tile) {
    if (id != MOLECULE_ID_INVALID) {
      molecule_ids.push_back(id);
    }
  }
}


void Grid::dump() const {
  // dumping just occupied locations and base info for now
  cout << "Grid: num_tiles: " << num_tiles << ", num_occupied: " << num_occupied << "\n";
  for (uint i = 0; i < molecules_per_tile.size(); i++) {
    molecule_id_t id = molecules_per_tile[i];
    if (id != MOLECULE_ID_INVALID) {
      cout << "[" << i << "]" << id << "\n";
    }
  }
}


/***************************************************************************
init_edge_transform
  In: e: pointer to an edge
      edgenum: integer telling which edge (0-2) of the "forward" face we are
  Out: No return value.  Coordinate transform in edge struct is set.
  Note: Don't call this on a non-shared edge.
***************************************************************************/
void Edge::reinit_edge_constants(const Partition& p) {
  assert(is_shared_edge());

  if (!is_initialized()) {
    return;
  }

  Wall wf = p.get_wall(forward_index);
  Wall wb = p.get_wall(backward_index);

#ifdef DEBUG_EDGE_INITIALIZATION
  std::cout << "Edge initialization, edgenum: " << edge_num_used_for_init << "\n";
  wf.dump(p, "", true);
  wb.dump(p, "", true);
#endif

  edge_index_t i = edge_num_used_for_init;
  assert(i < VERTICES_IN_TRIANGLE);
  edge_index_t j = i + 1;
  if (j == VERTICES_IN_TRIANGLE)
    j = 0;

  const Vec3& wf_vert_0 = p.get_geometry_vertex(wf.vertex_indices[0]);
  const Vec3& wf_vert_i = p.get_geometry_vertex(wf.vertex_indices[i]);
  const Vec3& wf_vert_j = p.get_geometry_vertex(wf.vertex_indices[j]);

  const Vec3& wb_vert_0 = p.get_geometry_vertex(wb.vertex_indices[0]);


  /* Intermediate basis from the perspective of the forward frame */
  Vec3 diff_i_0 = wf_vert_i - wf_vert_0;

  Vec2 O_f;
  O_f.u = dot(diff_i_0, wf.unit_u); // should be called after wall init
  O_f.v = dot(diff_i_0, wf.unit_v); /* Origin */

  Vec3 diff_j_0 = wf_vert_j - wf_vert_0;
  Vec2 temp_ff;
  temp_ff.u = dot(diff_j_0, wf.unit_u) - O_f.u;
  temp_ff.v = dot(diff_j_0, wf.unit_v) - O_f.v; /* Far side of e */

  assert(!cmp_eq(len2_squared(temp_ff), 0.0, EPS));
  float_t d_f = 1.0 / len2(temp_ff);

  Vec2 ehat_f, fhat_f;
  ehat_f = temp_ff * d_f; /* ehat along edge */
  fhat_f.u = -ehat_f.v;
  fhat_f.v = ehat_f.u; /* fhat 90 degrees CCW */

  /* Intermediate basis from the perspective of the backward frame */
  Vec3 diff_i_b0 = wf_vert_i - wb_vert_0;
  Vec2 O_b;
  O_b.u = dot(diff_i_b0, wb.unit_u);
  O_b.v = dot(diff_i_b0, wb.unit_v); /* Origin */

  Vec3 diff_j_b0 = wf_vert_j - wb_vert_0;
  Vec2 temp_fb;
  temp_fb.u = dot(diff_j_b0, wb.unit_u) - O_b.u;
  temp_fb.v = dot(diff_j_b0, wb.unit_v) - O_b.v; /* Far side of e */

  assert(!cmp_eq(len2_squared(temp_fb), 0.0, EPS));
  float_t d_b = 1.0 / len2(temp_fb);

  Vec2 ehat_b, fhat_b;

  ehat_b = temp_fb * d_b; /* ehat along edge */
  fhat_b.u = -ehat_b.v;
  fhat_b.v = ehat_b.u; /* fhat 90 degrees CCW */

  /* Calculate transformation matrix */

  float_t mtx[2][2];
  mtx[0][0] = ehat_f.u * ehat_b.u + fhat_f.u * fhat_b.u;
  mtx[0][1] = ehat_f.v * ehat_b.u + fhat_f.v * fhat_b.u;
  mtx[1][0] = ehat_f.u * ehat_b.v + fhat_f.u * fhat_b.v;
  mtx[1][1] = ehat_f.v * ehat_b.v + fhat_f.v * fhat_b.v;

  /* Calculate translation vector */

  Vec2 q;
  q = O_b;

  q.u -= mtx[0][0] * O_f.u + mtx[0][1] * O_f.v;
  q.v -= mtx[1][0] * O_f.u + mtx[1][1] * O_f.v;

  /* Store the results */
  cos_theta = mtx[0][0];
  sin_theta = mtx[0][1];
  translate = q;

#ifdef DEBUG_EDGE_INITIALIZATION
  dump();
#endif
}

void Edge::debug_check_values_are_uptodate(const Partition& p) {
  if (!is_initialized()) {
    return;
  }
#ifndef NDEBUG
  Vec2 orig_translate = translate;
#endif
  float_t orig_cos_theta = cos_theta;
  float_t orig_sin_theta = sin_theta;
#ifdef DEBUG_EDGE_INITIALIZATION
  dump();
#endif
  reinit_edge_constants(p);
#ifdef DEBUG_EDGE_INITIALIZATION
  dump();
#endif
  assert(cmp_eq(orig_translate, translate));
  assert(cmp_eq(orig_cos_theta, cos_theta));
  assert(cmp_eq(orig_sin_theta, sin_theta));
}


void Edge::dump(const std::string ind) const {
  cout << ind <<
      "Edge: translate: " << translate <<
      ", cos_theta: " << cos_theta <<
      ", sin_theta: " << sin_theta <<
      ", edge_num_used_for_init: " << edge_num_used_for_init << "\n";
}


void Wall::precompute_wall_constants(const Partition& p) {

  const Vec3& v0 = p.get_geometry_vertex(vertex_indices[0]);
  const Vec3& v1 = p.get_geometry_vertex(vertex_indices[1]);
  const Vec3& v2 = p.get_geometry_vertex(vertex_indices[2]);

  Vec3 vA, vB, vX;
  vA = v1 - v0;
  vB = v2 - v0;
  vX = cross(vA, vB);

  area = 0.5 * len3(vX);

  if (!distinguishable(area, 0, EPS)) {
    /* this is a degenerate polygon.
    * perform initialization and quit. */
    normal = Vec3(0);
    unit_u = Vec3(0);
    unit_v = Vec3(0);
    uv_vert1_u = 0;
    uv_vert2 = Vec2(0);
    distance_to_origin = 0;
    return;
  }

  Vec3 f1 = v1 - v0;
  float_t f1_len_squared = len3_squared(f1);
  assert(f1_len_squared != 0);
  float_t inv_f1_len = 1 / sqrt_f(f1_len_squared);

  unit_u = f1 * Vec3(inv_f1_len);

  Vec3 f2 = v2 - v0;
  normal = cross(unit_u, f2);

  float_t norm_len_squared = len3_squared(normal);
  assert(norm_len_squared != 0);
  float_t inv_norm_len = 1 / sqrt_f(norm_len_squared);
  normal = normal * Vec3(inv_norm_len);

  unit_v = cross(normal, unit_u);

  distance_to_origin =  dot(v0, normal);

  uv_vert1_u = dot(f1, unit_u);
  uv_vert2.u = dot(f2, unit_u);
  uv_vert2.v = dot(f2, unit_v);

  wall_constants_precomputed = true;
}

void Wall::reinit_edge_constants(const Partition& p) {
  // edges, uses info set by init_tri_wall
  for (uint edge_index = 0; edge_index < EDGES_IN_TRIANGLE; edge_index++) {
    if (edges[edge_index].is_shared_edge()) {
      edges[edge_index].reinit_edge_constants(p);
    }
  }
}


void Wall::dump(const Partition& p, const std::string ind, const bool for_diff) const {
  if (for_diff) {
    cout << "wall[side: " << side << "]: ";
    for (uint i = 0; i < VERTICES_IN_TRIANGLE; i++) {
      vertex_index_t vertex_index = vertex_indices[i];
      Vec3 pos = p.get_geometry_vertex(vertex_index);
      cout << pos;
      if (i != VERTICES_IN_TRIANGLE - 1) {
        cout << ", ";
      }
    }
    cout << "\n";
  }
  else {
    cout << ind << "id: " << id << ", side: " << side << ", object_id: " << object_id << "\n";

    for (uint i = 0; i < VERTICES_IN_TRIANGLE; i++) {
      vertex_index_t vertex_index = vertex_indices[i];
      Vec3 pos = p.get_geometry_vertex(vertex_index);
      cout << ind << "vertex_index: " << vertex_index << ":" << pos << "\n";
    }

    for (uint i = 0; i < EDGES_IN_TRIANGLE; i++) {
      cout << ind << "edges:\n";
      edges[i].dump("        ");
    }

    cout << ind;
    for (uint i = 0; i < EDGES_IN_TRIANGLE; i++) {
      cout << "nb_walls[" << i << "]: " << nb_walls[i] << ", ";
    }
    cout << "\n";

    cout << ind << "regions: ";
    for (region_index_t i: regions) {
      cout << i << ", ";
    }
    cout << "\n";

    cout << ind
      << "uv_vert1_u: " << uv_vert1_u
      << ", uv_vert2: " << uv_vert2
      << ", area: " << area
      << ", normal: " << normal << "\n";

    cout << ind
      << "unit_u: " << unit_u
      << ", unit_v: " << unit_v
      << ", distance_to_origin: " << distance_to_origin
      << "\n";
  }
}


namespace Geometry {

// this is the entry point called from Partition class
void update_moved_walls(
    Partition& p,
    const std::vector<VertexMoveInfo>& scheduled_vertex_moves,
    // we can compute all the information already from scheduled_vertex_moves,
    // but the keys of the map walls_with_their_moves are the walls that we need to update
    const WallsWithTheirMovesMap& walls_with_their_moves
) {

  // move all vertices
  for (const VertexMoveInfo& move_info: scheduled_vertex_moves) {
    Vec3& vertex_ref = p.get_geometry_vertex(move_info.vertex_index);
    vertex_ref = vertex_ref + move_info.translation_vec;
    if (! p.in_this_partition(vertex_ref) ) {
      mcell_log("Error: Crossing partitions is not supported yet.\n");
      exit(1);
    }
  }

  // update walls
  // FIXME: this should be placed in geometry.cpp
  for (auto it: walls_with_their_moves) {
    wall_index_t wall_index = it.first;
    Wall& w = p.get_wall(wall_index);

    // first we need to update all wall constants
    w.precompute_wall_constants(p);

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
    w.reinit_edge_constants(p);
  }
}

} /* namespace Geometry */

} /* namespace MCell */
