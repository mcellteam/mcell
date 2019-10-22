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

void Grid::initialize(const Partition& p, const Wall& w) {

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


void GeometryObject::dump(const Partition& p, const std::string ind) const {
  cout << ind << "geometry_object_t: id:" << id << ", name:" << name << "\n";
  for (wall_index_t i: wall_indices) {
    cout << ind << "  " << i << ": ";
    p.get_wall(i).dump(p, ind + "  ");
  }
}


// should belong to Edge?
static void init_edge_transform(Partition& p, Edge &e, int edgenum) {

  // not sure what to do if these asserts do not hold
  assert(e.forward_index != WALL_INDEX_INVALID);
  assert(e.backward_index != WALL_INDEX_INVALID);
  Wall wf = p.get_wall(e.forward_index);
  Wall wb = p.get_wall(e.backward_index);

  uint i = edgenum;
  assert(i < VERTICES_IN_TRIANGLE);
  uint j = i + 1;
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

  vec3_t diff_j_b0 = wf_vert_i - wb_vert_0;
  vec2_t temp_fb;
  temp_fb.u = dot(diff_j_b0, wb.unit_u) - O_b.u;
  temp_fb.v = dot(diff_j_b0, wb.unit_v) - O_b.v; /* Far side of e */

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
  e.cos_theta = mtx[0][0];
  e.sin_theta = mtx[0][1];
  e.translate = q;
}



static void init_tri_wall(Partition& p, Wall& w) {

  const vec3_t& v0 = p.get_geometry_vertex(w.vertex_indices[0]);
  const vec3_t& v1 = p.get_geometry_vertex(w.vertex_indices[1]);
  const vec3_t& v2 = p.get_geometry_vertex(w.vertex_indices[2]);

  vec3_t vA, vB, vX;
  vA = v1 - v0;
  vB = v2 - v0;
  vX = cross(vA, vB);

  w.area = 0.5 * len3(vX); // vect_length(&vX);

  if (!distinguishable(w.area, 0, EPS_C)) {
    /* this is a degenerate polygon.
    * perform initialization and quit. */
    w.normal = vec3_t(0);
    w.unit_u = vec3_t(0);
    w.unit_v = vec3_t(0);
    w.uv_vert1_u = 0;
    w.uv_vert2 = vec2_t(0);
    w.distance_to_origin = 0;
    return;
  }

  // FIXME: simplify
  float_t f, fx, fy, fz;
  fx = (v1.x - v0.x);
  fy = (v1.y - v0.y);
  fz = (v1.z - v0.z);
  f = 1 / sqrt(fx * fx + fy * fy + fz * fz);

  w.unit_u.x = fx * f;
  w.unit_u.y = fy * f;
  w.unit_u.z = fz * f;

  fx = (v2.x - v0.x);
  fy = (v2.y - v0.y);
  fz = (v2.z - v0.z);

  w.normal.x = w.unit_u.y * fz - w.unit_u.z * fy;
  w.normal.y = w.unit_u.z * fx - w.unit_u.x * fz;
  w.normal.z = w.unit_u.x * fy - w.unit_u.y * fx;
  f = 1 / sqrt(w.normal.x * w.normal.x + w.normal.y * w.normal.y +
               w.normal.z * w.normal.z);
  w.normal.x *= f;
  w.normal.y *= f;
  w.normal.z *= f;
  w.unit_v.x = w.normal.y * w.unit_u.z - w.normal.z * w.unit_u.y;
  w.unit_v.y = w.normal.z * w.unit_u.x - w.normal.x * w.unit_u.z;
  w.unit_v.z = w.normal.x * w.unit_u.y - w.normal.y * w.unit_u.x;
  w.distance_to_origin = v0.x * w.normal.x + v0.y * w.normal.y + v0.z * w.normal.z;

  w.uv_vert1_u = (v1.x - v0.x) * w.unit_u.x +
                  (v1.y - v0.y) * w.unit_u.y +
                  (v1.z - v0.z) * w.unit_u.z;
  w.uv_vert2.u = (v2.x - v0.x) * w.unit_u.x +
                  (v2.y - v0.y) * w.unit_u.y +
                  (v2.z - v0.z) * w.unit_u.z;
  w.uv_vert2.v = (v2.x - v0.x) * w.unit_v.x +
                  (v2.y - v0.y) * w.unit_v.y +
                  (v2.z - v0.z) * w.unit_v.z;
}


void Wall::update_after_vertex_change(Partition& p) {

  // the only thing that we know that is correct now is the
  // position of the vertices, we need to recompute everything else

  // area
  // normal
  // unit_u
  // unit_v
  // distance_to_origin
  init_tri_wall(p, *this);

  // edges, uses info set by init_tri_wall
  for (uint edge_index = 0; edge_index < EDGES_IN_TRIANGLE; edge_index++) {
    assert(edges[edge_index].forward_index == index); // ???
    init_edge_transform(p, edges[edge_index], edge_index /*???*/);
  }

  // walls must not have any surface molecules for now!
  // TODO: do we need a better check?
#ifndef NDEBUG
  for (tile_index_t tile_index = 0; tile_index < grid.num_tiles; tile_index++) {
    assert(grid.get_molecule_on_tile(tile_index) == MOLECULE_INDEX_INVALID);
  }
#endif

  // reinitialize grid
  grid.initialize(p, *this);
}


void Wall::dump(const Partition& p, const std::string ind) const {
  cout << "id: " << id << ", side: " << side << ", object_id: " << object_id;

  for (uint i = 0; i < VERTICES_IN_TRIANGLE; i++) {
    vertex_index_t vertex_index = vertex_indices[i];
    vec3_t pos = p.get_geometry_vertex(vertex_index);
    cout << ", vert_index: " << vertex_index << ":" << pos;
  }

  cout
    << ", normal " << normal
    << ", unit_u " << unit_u
    << ", unit_v " << unit_v
    << "\n";
}

} /* namespace mcell */
