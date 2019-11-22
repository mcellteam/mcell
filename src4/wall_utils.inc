/******************************************************************************
 *
 * Copyright (C) 2006-2017 by
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

#ifndef SRC4_WALL_UTILS_INC_
#define SRC4_WALL_UTILS_INC_

/**
 * This file is directly included into diffuse_react_event.cpp.
 * The reason why this is not a standard .cpp + .h file is to gove the compiler
 * the opportunity to inline these functions into methods of diffuse&react event.
 */
#include <vector>

#include "diffuse_react_event.h"
#include "defines.h"
#include "world.h"
#include "partition.h"
#include "geometry.h"
#include "debug_config.h"

namespace MCell {

namespace WallUtil {

/***********************************************************************
walls_share_full_edge:
  In: two walls
  Out: 1 if the walls share a full edge, 0 - otherwise.
       Here by "full" we mean that the shared edge has two endpoints
       that are the vertices of both walls w1 and w2.
************************************************************************/
static bool walls_share_full_edge(const Partition& p, const Wall& w1, const Wall& w2) {
  uint i, k;
  uint count = 0;

  /* count number of shared vertices between two walls */
  for (i = 0; i < VERTICES_IN_TRIANGLE; i++) {
    for (k = 0; k < VERTICES_IN_TRIANGLE; k++) {

      if (!distinguishable_vec3( p.get_wall_vertex(w1, i), p.get_wall_vertex(w2, k), EPS) ) {
        count++;
      }
    }
  }

  return count == 2;
}

/*************************************************************************
find_nbr_walls_shared_one_vertex:
   In: the origin wall
       array with information about which vertices of the origin wall
          are shared with neighbor wall (they are indices in the
          global "world->walls_using_vertex" array).
   Out: linked list of the neighbor walls that have only one common
        vertex with the origin wall (not edge-to-edge walls, but
        vertex-to-vertex walls).
   Note: the "origin" wall is not included in the list
**************************************************************************/
static void find_nbr_walls_shared_one_vertex(
    const Partition& p,
    const Wall& origin_wall,
    vertex_index_t shared_verts[VERTICES_IN_TRIANGLE],
    wall_indices_t& neighboring_walls
)
{
  for (uint i = 0; i < VERTICES_IN_TRIANGLE; i++) {
    if (shared_verts[i] != VERTEX_INDEX_INVALID) {

      const std::vector<wall_index_t>& wall_indices = p.get_walls_using_vertex(shared_verts[i]);
      for (wall_index_t wi: wall_indices) {
        const Wall& w = p.get_wall(wi);

        if (w.id == origin_wall.id) {
          // we do not care about current wall
          continue;
        }

        if (!walls_share_full_edge(p, origin_wall, w)) {
          neighboring_walls.push_back(w.id);
        }
      }
    }
  }
}



/***************************************************************************
wall_in_box:
  In: array of pointers to vertices for wall (should be 3)
      normal vector for wall
      distance from wall to origin (point normal form)
      first corner of bounding box
      opposite corner of bounding box
  Out: nonzero if the wall intersects the box.  0 otherwise.
***************************************************************************/
// TODO: cleanup
static int wall_in_box(
    const Partition& p,
    const Wall& w,
    const vec3_t& llf, const vec3_t& urb) {

  const vec3_t* vert[VERTICES_IN_TRIANGLE] = {
      &p.get_wall_vertex(w, 0),
      &p.get_wall_vertex(w, 1),
      &p.get_wall_vertex(w, 2)
  };
  const vec3_t& normal = w.normal;
  float_t d = w.distance_to_origin;

  int temp;
  uint i, j, k;
  const vec3_t *v1, *v2;
  vec3_t n, u, v;
  vec3_t ba, bb, c;
  float_t r, a1, a2, a3, a4, cu, cv;
  float_t vu_[VERTICES_IN_TRIANGLE * 2]; /* Assume wall has 3 vertices */
  float_t *vv_;
  float_t d_box[8];
  int n_opposite;

  /* Lookup table for vertex-edge mapping for a cube */
  int which_x1[12] = { 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1 };
  int which_y1[12] = { 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 0 };
  int which_z1[12] = { 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0 };
  int which_x2[12] = { 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0 };
  int which_y2[12] = { 0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1 };
  int which_z2[12] = { 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0 };

  int edge1_vt[12] = { 0, 1, 3, 2, 6, 7, 5, 4, 0, 1, 3, 4 };
  int edge2_vt[12] = { 1, 3, 2, 6, 7, 5, 4, 0, 2, 5, 7, 2 };

  /* Check if any vertex of the wall is in the box. */
  for (i = 0; i < VERTICES_IN_TRIANGLE; i++) {
    v2 = vert[i];
    if (v2->x >= llf.x && v2->x <= urb.x && v2->y >= llf.y && v2->y <= urb.y &&
        v2->z >= llf.z && v2->z <= urb.z)
      return 1;
  }

  /* Check if any wall edge intersects any face of the box */
  for (i = 0; i < VERTICES_IN_TRIANGLE; i++) {
    v2 = vert[i];
    v1 = (i == 0) ? vert[VERTICES_IN_TRIANGLE - 1] : vert[i - 1];

    /* x-faces */
    if ((v1->x <= llf.x && llf.x < v2->x) ||
        (v1->x > llf.x && llf.x >= v2->x)) {
      r = (llf.x - v1->x) / (v2->x - v1->x);
      a3 = v1->y + r * (v2->y - v1->y);
      a4 = v1->z + r * (v2->z - v1->z);
      if (llf.y <= a3 && a3 <= urb.y && llf.z <= a4 && a4 <= urb.z)
        return 2;
    }
    if ((v1->x <= urb.x && urb.x < v2->x) ||
        (v1->x > urb.x && urb.x >= v2->x)) {
      r = (urb.x - v1->x) / (v2->x - v1->x);
      a3 = v1->y + r * (v2->y - v1->y);
      a4 = v1->z + r * (v2->z - v1->z);
      if (llf.y <= a3 && a3 <= urb.y && llf.z <= a4 && a4 <= urb.z)
        return 3;
    }

    /* y-faces */
    if ((v1->y <= llf.y && llf.y < v2->y) ||
        (v1->y > llf.y && llf.y >= v2->y)) {
      r = (llf.y - v1->y) / (v2->y - v1->y);
      a3 = v1->x + r * (v2->x - v1->x);
      a4 = v1->z + r * (v2->z - v1->z);
      if (llf.x <= a3 && a3 <= urb.x && llf.z <= a4 && a4 <= urb.z)
        return 4;
    }
    if ((v1->y <= urb.y && urb.y < v2->y) ||
        (v1->y > urb.y && urb.y >= v2->y)) {
      r = (urb.y - v1->y) / (v2->y - v1->y);
      a3 = v1->x + r * (v2->x - v1->x);
      a4 = v1->z + r * (v2->z - v1->z);
      if (llf.x <= a3 && a3 <= urb.x && llf.z <= a4 && a4 <= urb.z)
        return 5;
    }

    /* z-faces */
    if ((v1->z <= llf.z && llf.z < v2->z) ||
        (v1->z > llf.z && llf.z >= v2->z)) {
      r = (llf.z - v1->z) / (v2->z - v1->z);
      a3 = v1->y + r * (v2->y - v1->y);
      a4 = v1->x + r * (v2->x - v1->x);
      if (llf.y <= a3 && a3 <= urb.y && llf.x <= a4 && a4 <= urb.x)
        return 6;
    }
    if ((v1->z <= urb.z && urb.z < v2->z) ||
        (v1->z > urb.z && urb.z >= v2->z)) {
      r = (urb.z - v1->z) / (v2->z - v1->z);
      a3 = v1->y + r * (v2->y - v1->y);
      a4 = v1->x + r * (v2->x - v1->x);
      if (llf.y <= a3 && a3 <= urb.y && llf.x <= a4 && a4 <= urb.x)
        return 7;
    }
  }

  /* Check if any box edge intersects the wall */

  n_opposite = 0;
  vv_ = &(vu_[VERTICES_IN_TRIANGLE]);

  /* Wall coordinate system n,u,v */
  n.x = normal.x;
  n.y = normal.y;
  n.z = normal.z;
  u.x = vert[1]->x - vert[0]->x;
  u.y = vert[1]->y - vert[0]->y;
  u.z = vert[1]->z - vert[0]->z;
  r = 1 / sqrt(u.x * u.x + u.y * u.y + u.z * u.z);
  u.x *= r;
  u.y *= r;
  u.z *= r;
  v.x = n.y * u.z - n.z * u.y;
  v.y = -(n.x * u.z - n.z * u.x);
  v.z = n.x * u.y - n.y * u.x;
  for (j = 0; j < VERTICES_IN_TRIANGLE; j++) {
    vu_[j] = vert[j]->x * u.x + vert[j]->y * u.y + vert[j]->z * u.z;
    vv_[j] = vert[j]->x * v.x + vert[j]->y * v.y + vert[j]->z * v.z;
  }

  /* Test every edge. */
  bb.x = llf.x;
  bb.y = llf.y;
  bb.z = llf.z;
  d_box[0] = bb.x * n.x + bb.y * n.y + bb.z * n.z;
  for (i = 0; i < 12; i++) {
    if (i < 7) /* Visiting new vertices in order */
    {
      ba.x = bb.x;
      ba.y = bb.y;
      ba.z = bb.z;
      bb.x = (which_x2[i]) ? urb.x : llf.x;
      bb.y = (which_y2[i]) ? urb.y : llf.y;
      bb.z = (which_z2[i]) ? urb.z : llf.z;
      a2 = d_box[edge2_vt[i]] = bb.x * n.x + bb.y * n.y + bb.z * n.z;
      a1 = d_box[edge1_vt[i]];

      if ((a1 - d < 0 && a2 - d < 0) || (a1 - d > 0 && a2 - d > 0))
        continue;
      else
        n_opposite++;
    } else /* Revisiting old vertices out of order */
    {
      /*      if (!n_opposite) return 0; */
      a1 = d_box[edge1_vt[i]];
      a2 = d_box[edge2_vt[i]];

      if ((a1 - d < 0 && a2 - d < 0) || (a1 - d > 0 && a2 - d > 0))
        continue;

      n_opposite++;
      ba.x = (which_x1[i]) ? urb.x : llf.x;
      ba.y = (which_y1[i]) ? urb.y : llf.y;
      ba.z = (which_z1[i]) ? urb.z : llf.z;
      bb.x = (which_x2[i]) ? urb.x : llf.x;
      bb.y = (which_y2[i]) ? urb.y : llf.y;
      bb.z = (which_z2[i]) ? urb.z : llf.z;
    }
    /* Now ba,bb = box edge endpoints ; a1,a2 = distances along wall normal */
    r = (d - a1) / (a2 - a1);
    c.x = ba.x + r * (bb.x - ba.x);
    c.y = ba.y + r * (bb.y - ba.y);
    c.z = ba.z + r * (bb.z - ba.z);
    cu = c.x * u.x + c.y * u.y + c.z * u.z;
    cv = c.x * v.x + c.y * v.y + c.z * v.z;
    /* Test for internal intersection point in wall coordinate space */
    temp = 0;
    for (j = 0; j < VERTICES_IN_TRIANGLE; j++) {
      k = (j == 0) ? VERTICES_IN_TRIANGLE - 1 : j - 1;
      if ((vu_[k] < cu && cu <= vu_[j]) || (vu_[k] >= cu && cu > vu_[j])) {
        r = (cu - vu_[k]) / (vu_[j] - vu_[k]);
        if ((vv_[k] + r * (vv_[j] - vv_[k])) > cv)
          temp++;
      }
    }
    if (temp & 1)
      return 8 + i;
  }

  return 0;
}


} // namespace wall_util

} // namespace mcell

#endif // SRC4_WALL_UTILS_INC_
