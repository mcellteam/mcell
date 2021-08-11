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

using namespace std;

namespace MCell {

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

  num_tiles_along_axis = (int)ceil_p(sqrt_p(w.area));
  if (num_tiles_along_axis < 1) {
    num_tiles_along_axis = 1;
  }

  num_tiles = num_tiles_along_axis * num_tiles_along_axis;

  molecules_per_tile.resize(num_tiles, MOLECULE_ID_INVALID);

  strip_width_rcp = 1 / (w.uv_vert2.v / ((pos_t)num_tiles_along_axis));
  vert2_slope = w.uv_vert2.u / w.uv_vert2.v;
  fullslope = w.uv_vert1_u / w.uv_vert2.v;

  binding_factor = ((pos_t)num_tiles) / w.area;

  const Vec3& vert0_tmp = p.get_wall_vertex(w, 0);

  vert0.u = dot(vert0_tmp, w.unit_u);
  vert0.v = dot(vert0_tmp, w.unit_v);

  assert(w.index != WALL_INDEX_INVALID);
  wall_index = w.index;
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


static void report_malformed_geom_object_and_exit(const Partition& p, const Wall& wf, const Wall& wb) {
  const GeometryObject& o = p.get_geometry_object(wb.object_index);
  errs() << "Detected malformed geometry object '" << o.name << "'. " <<
      "A possible cause is that two vertices share the same location or a wall became very thin (like a line).\n" <<
      "Error detected for walls with side indices " << wf.side << " and " << wb.side << ".\n" <<
      "Terminating because this issue would lead to simulation errors.\n";
  exit(1);
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

  if (cmp_eq(len2_squared(temp_ff), (pos_t)0, POS_EPS)) {
    report_malformed_geom_object_and_exit(p, wf, wb);
  }

  pos_t d_f = 1 / len2(temp_ff);

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

  if (cmp_eq(len2_squared(temp_fb), (pos_t)0, POS_EPS)) {
    report_malformed_geom_object_and_exit(p, wf, wb);
  }

  pos_t d_b = 1 / len2(temp_fb);

  Vec2 ehat_b, fhat_b;

  ehat_b = temp_fb * d_b; /* ehat along edge */
  fhat_b.u = -ehat_b.v;
  fhat_b.v = ehat_b.u; /* fhat 90 degrees CCW */

  /* Calculate transformation matrix */

  pos_t mtx[2][2];
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
  pos_t orig_cos_theta = cos_theta;
  pos_t orig_sin_theta = sin_theta;
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
  Vec2 translate_rounded;
  translate_rounded.x = (translate.x < POS_EPS) ? 0 : translate.x;
  translate_rounded.y = (translate.y < POS_EPS) ? 0 : translate.y;

  cout << ind <<
      "Edge: "
      "forward_index: " << forward_index <<
      ", backward_index: " << backward_index <<
      ", translate: " << translate_rounded <<
      ", cos_theta: " << ((cos_theta < POS_EPS) ? 0 : cos_theta) <<
      ", sin_theta: " << ((sin_theta < POS_EPS) ? 0 : sin_theta) <<
      ", edge_num_used_for_init: " << edge_num_used_for_init << "\n";
}


void Wall::initalize_wall_constants(const Partition& p) {

  const Vec3* v0;
  const Vec3* v1;
  const Vec3* v2;
  if (exists_in_partition()) {
    v0 = &p.get_geometry_vertex(vertex_indices[0]);
    v1 = &p.get_geometry_vertex(vertex_indices[1]);
    v2 = &p.get_geometry_vertex(vertex_indices[2]);
  }
  else {
    WallWithVertices* wall_w_vertices = static_cast<WallWithVertices*>(this);
    assert(wall_w_vertices != nullptr);
    v0 = &wall_w_vertices->vertices[0];
    v1 = &wall_w_vertices->vertices[1];
    v2 = &wall_w_vertices->vertices[2];
  }

  Vec3 vA, vB, vX;
  vA = *v1 - *v0;
  vB = *v2 - *v0;
  vX = cross(vA, vB);

  area = 0.5 * len3(vX);

  if (!distinguishable_f(area, 0, EPS)) {
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

  Vec3 f1 = *v1 - *v0;
  pos_t f1_len_squared = len3_squared(f1);
  assert(f1_len_squared != 0);
  pos_t inv_f1_len = 1 / sqrt_p(f1_len_squared);

  unit_u = f1 * Vec3(inv_f1_len);

  Vec3 f2 = *v2 - *v0;
  normal = cross(unit_u, f2);

  pos_t norm_len_squared = len3_squared(normal);
  assert(norm_len_squared != 0);
  pos_t inv_norm_len = 1 / sqrt_p(norm_len_squared);
  normal = normal * Vec3(inv_norm_len);

  unit_v = cross(normal, unit_u);

  distance_to_origin =  dot(*v0, normal);

  uv_vert1_u = dot(f1, unit_u);
  uv_vert2.u = dot(f2, unit_u);
  uv_vert2.v = dot(f2, unit_v);

  wall_constants_initialized = true;
}

void Wall::initialize_edge_constants(const Partition& p) {
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
    cout << ind <<
        "id: " << id << ", index: " << index << ", side: " << side <<
        ", object_id: " << object_id << ", object_index: " << object_index << "\n";

    for (uint i = 0; i < VERTICES_IN_TRIANGLE; i++) {
      vertex_index_t vertex_index = vertex_indices[i];
      Vec3 pos = p.get_geometry_vertex(vertex_index);
      cout << ind << "vertex_index: " << vertex_index << ":" << pos << "\n";
    }

    cout << ind << "edges:\n";
    for (uint i = 0; i < EDGES_IN_TRIANGLE; i++) {
      cout << ind << i << ":\n";
      edges[i].dump(ind + "  ");
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


/***************************************************************************
create_region_bbox:
  In: a region
  Out: LLF and URB corners of a bounding box around the region
***************************************************************************/
void Region::compute_bounding_box(const Partition& p, Vec3& llf, Vec3& urb) {
  bool found_first_wall = false;
  for (auto it: walls_and_edges) {
    const Wall& w = p.get_wall(it.first);
    const Vec3* verts[VERTICES_IN_TRIANGLE];
    verts[0] = &p.get_wall_vertex(w, 0);
    verts[1] = &p.get_wall_vertex(w, 1);
    verts[2] = &p.get_wall_vertex(w, 2);

    if (!found_first_wall) {
      llf = *verts[0];
      urb = *verts[0];
      found_first_wall = true;
    }

    for (uint i = 0; i < VERTICES_IN_TRIANGLE; i++) {
      const Vec3* v = verts[i];
      if (llf.x > v->x) {
        llf.x = v->x;
      }
      else if (urb.x < v->x) {
        urb.x = v->x;
      }

      if (llf.y > v->y) {
        llf.y = v->y;
      }
      else if (urb.y < v->y) {
        urb.y = v->y;
      }

      if (llf.z > v->z) {
        llf.z = v->z;
      }
      else if (urb.z < v->z) {
        urb.z = v->z;
      }
    }
  }
}

} /* namespace MCell */
