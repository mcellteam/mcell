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

//extern "C" {
#include "rng.h" // MCell 3
#include "isaac64.h"
#include "mcell_structs.h"
#include "logging.h"
//}

#include <iostream>

#include "partition.h"
#include "geometry.h"

#include "geometry_utils.inc" // uses get_wall_bounding_box, maybe not include this file

using namespace std;

namespace MCell {

// may be also used fore reinitialization
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

  const vec3_t& vert0_tmp = p.get_wall_vertex(w, 0);

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

void GeometryObject::dump(const Partition& p, const std::string ind) const {
  cout << ind << "geometry_object_t: id:" << id << ", name:" << name << "\n";
  for (wall_index_t i: wall_indices) {
    cout << ind << "  " << i << ": \n";
    p.get_wall(i).dump(p, ind + "    ");
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

  if (!is_initialized()) {
    return;
  }

  // not sure what to do if these asserts do not hold
  assert(forward_index != WALL_INDEX_INVALID);
  assert(backward_index != WALL_INDEX_INVALID);
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

  const vec3_t& wf_vert_0 = p.get_geometry_vertex(wf.vertex_indices[0]);
  const vec3_t& wf_vert_i = p.get_geometry_vertex(wf.vertex_indices[i]);
  const vec3_t& wf_vert_j = p.get_geometry_vertex(wf.vertex_indices[j]);

  const vec3_t& wb_vert_0 = p.get_geometry_vertex(wb.vertex_indices[0]);


  /* Intermediate basis from the perspective of the forward frame */
  vec3_t diff_i_0 = wf_vert_i - wf_vert_0;

  vec2_t O_f;
  O_f.u = dot(diff_i_0, wf.unit_u); // should be called after wall init
  O_f.v = dot(diff_i_0, wf.unit_v); /* Origin */

  vec3_t diff_j_0 = wf_vert_j - wf_vert_0;
  vec2_t temp_ff;
  temp_ff.u = dot(diff_j_0, wf.unit_u) - O_f.u;
  temp_ff.v = dot(diff_j_0, wf.unit_v) - O_f.v; /* Far side of e */

  assert(!cmp_eq(len2_squared(temp_ff), 0.0, EPS));
  float_t d_f = 1.0 / len2(temp_ff);

  vec2_t ehat_f, fhat_f;
  ehat_f = temp_ff * d_f; /* ehat along edge */
  fhat_f.u = -ehat_f.v;
  fhat_f.v = ehat_f.u; /* fhat 90 degrees CCW */

  /* Intermediate basis from the perspective of the backward frame */
  vec3_t diff_i_b0 = wf_vert_i - wb_vert_0;
  vec2_t O_b;
  O_b.u = dot(diff_i_b0, wb.unit_u);
  O_b.v = dot(diff_i_b0, wb.unit_v); /* Origin */

  vec3_t diff_j_b0 = wf_vert_j - wb_vert_0;
  vec2_t temp_fb;
  temp_fb.u = dot(diff_j_b0, wb.unit_u) - O_b.u;
  temp_fb.v = dot(diff_j_b0, wb.unit_v) - O_b.v; /* Far side of e */

  assert(!cmp_eq(len2_squared(temp_fb), 0.0, EPS));
  float_t d_b = 1.0 / len2(temp_fb);

  vec2_t ehat_b, fhat_b;

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

  vec2_t q;
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
  if (edge_num_used_for_init != EDGE_INDEX_INVALID) {
    // not initialized
    return;
  }
  vec2_t orig_translate = translate;
  float_t orig_cos_theta = cos_theta;
  float_t orig_sin_theta = sin_theta;
  dump();
  reinit_edge_constants(p);
  dump();
  assert(cmp_eq(orig_translate, translate));
  assert(cmp_eq(orig_cos_theta, cos_theta));
  assert(cmp_eq(orig_sin_theta, sin_theta));
}

void Edge::dump() {
  cout <<
      "Edge: translate: " << translate <<
      ", cos_theta: " << cos_theta <<
      ", sin_theta: " << sin_theta <<
      ", edge_num_used_for_init: " << edge_num_used_for_init << "\n";
}


void Wall::precompute_wall_constants(const Partition& p) {

  const vec3_t& v0 = p.get_geometry_vertex(vertex_indices[0]);
  const vec3_t& v1 = p.get_geometry_vertex(vertex_indices[1]);
  const vec3_t& v2 = p.get_geometry_vertex(vertex_indices[2]);

  vec3_t vA, vB, vX;
  vA = v1 - v0;
  vB = v2 - v0;
  vX = cross(vA, vB);

  area = 0.5 * len3(vX); // vect_length(&vX);

  if (!distinguishable(area, 0, EPS_C)) {
    /* this is a degenerate polygon.
    * perform initialization and quit. */
    normal = vec3_t(0);
    unit_u = vec3_t(0);
    unit_v = vec3_t(0);
    uv_vert1_u = 0;
    uv_vert2 = vec2_t(0);
    distance_to_origin = 0;
    return;
  }

  // FIXME: simplify
  float_t f, fx, fy, fz;
  fx = (v1.x - v0.x);
  fy = (v1.y - v0.y);
  fz = (v1.z - v0.z);
  // TODO: assert
  f = 1 / sqrt(fx * fx + fy * fy + fz * fz);

  unit_u.x = fx * f;
  unit_u.y = fy * f;
  unit_u.z = fz * f;

  fx = (v2.x - v0.x);
  fy = (v2.y - v0.y);
  fz = (v2.z - v0.z);

  normal.x = unit_u.y * fz - unit_u.z * fy;
  normal.y = unit_u.z * fx - unit_u.x * fz;
  normal.z = unit_u.x * fy - unit_u.y * fx;
  // TODO: assert
  f = 1 / sqrt(normal.x * normal.x + normal.y * normal.y +
               normal.z * normal.z);
  normal.x *= f;
  normal.y *= f;
  normal.z *= f;
  unit_v.x = normal.y * unit_u.z - normal.z * unit_u.y;
  unit_v.y = normal.z * unit_u.x - normal.x * unit_u.z;
  unit_v.z = normal.x * unit_u.y - normal.y * unit_u.x;
  distance_to_origin = v0.x * normal.x + v0.y * normal.y + v0.z * normal.z;

  uv_vert1_u = (v1.x - v0.x) * unit_u.x +
                  (v1.y - v0.y) * unit_u.y +
                  (v1.z - v0.z) * unit_u.z;
  uv_vert2.u = (v2.x - v0.x) * unit_u.x +
                  (v2.y - v0.y) * unit_u.y +
                  (v2.z - v0.z) * unit_u.z;
  uv_vert2.v = (v2.x - v0.x) * unit_v.x +
                  (v2.y - v0.y) * unit_v.y +
                  (v2.z - v0.z) * unit_v.z;

  wall_constants_precomputed = true;
}

void Wall::reinit_edge_constants(const Partition& p) {
  // edges, uses info set by init_tri_wall
  for (uint edge_index = 0; edge_index < EDGES_IN_TRIANGLE; edge_index++) {
    edges[edge_index].reinit_edge_constants(p);
  }
}


void Wall::dump(const Partition& p, const std::string ind, const bool for_diff) const {
  if (for_diff) {
    cout << "wall[side: " << side << "]: ";
    for (uint i = 0; i < VERTICES_IN_TRIANGLE; i++) {
      vertex_index_t vertex_index = vertex_indices[i];
      vec3_t pos = p.get_geometry_vertex(vertex_index);
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
      vec3_t pos = p.get_geometry_vertex(vertex_index);
      cout << ind << "vertex_index: " << vertex_index << ":" << pos << "\n";
    }

    for (uint i = 0; i < EDGES_IN_TRIANGLE; i++) {
      cout << ind << "edges:\n";
      //TODO: edges[i].dump();
    }

    cout << ind;
    for (uint i = 0; i < EDGES_IN_TRIANGLE; i++) {
      cout << "nb_walls[" << i << "]: " << nb_walls[i] << ", ";
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

} /* namespace mcell */
