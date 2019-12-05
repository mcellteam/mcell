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

#ifndef SRC4_GEOMETRY_UTILS_INC_
#define SRC4_GEOMETRY_UTILS_INC_

/**
 * This file is directly included into diffuse_react_event.cpp.
 * The reason why this is not a standard .cpp + .h file is to gove the compiler
 * the opportunity to inline these functions into methods of diffuse&react event.
 */
#include <vector>

#include "logging.h"

#include "diffuse_react_event.h"
#include "defines.h"
#include "world.h"
#include "partition.h"
#include "geometry.h"
#include "debug_config.h"

#include "wall_utils.inc"

namespace MCell {

namespace GeometryUtil {

// several commonly used geometry utilities,
// rest of these utilities is in
static vec3_t uv2xyz(const vec2_t& a, const Wall& w, const vec3_t& wall_vert0) {
  return vec3_t(a.u) * w.unit_u + vec3_t(a.v) * w.unit_v + wall_vert0;
}


static vec2_t xyz2uv(const Partition& p, const vec3_t& a, const Wall& w) {
  vec2_t res;

  if (w.has_initialized_grid()) {
    const Grid& g = w.grid;

    res.u = a.x * w.unit_u.x + a.y * w.unit_u.y + a.z * w.unit_u.z -
           g.vert0.u;
    res.v = a.x * w.unit_v.x + a.y * w.unit_v.y + a.z * w.unit_v.z -
           g.vert0.v;
  }
  else {
    vec3_t pos = a - p.get_wall_vertex(w, 0);
    res.u = dot(pos, w.unit_u);
    res.v = dot(pos, w.unit_v);
  }
  return res;
}


/***************************************************************************
wall_bounding_box:
  In: a wall
      vector to store one corner of the bounding box for that wall
      vector to store the opposite corner
  Out: No return value.  The vectors are set to define the smallest box
       that contains the wall.
***************************************************************************/
static inline void get_wall_bounding_box(
    const vec3_t w_vert[VERTICES_IN_TRIANGLE],
    vec3_t& llf, vec3_t& urb
) {
  llf.x = urb.x = w_vert[0].x;
  llf.y = urb.y = w_vert[0].y;
  llf.z = urb.z = w_vert[0].z;

  if (w_vert[1].x < llf.x)
    llf.x = w_vert[1].x;
  else if (w_vert[1].x > urb.x)
    urb.x = w_vert[1].x;
  if (w_vert[2].x < llf.x)
    llf.x = w_vert[2].x;
  else if (w_vert[2].x > urb.x)
    urb.x = w_vert[2].x;

  if (w_vert[1].y < llf.y)
    llf.y = w_vert[1].y;
  else if (w_vert[1].y > urb.y)
    urb.y = w_vert[1].y;
  if (w_vert[2].y < llf.y)
    llf.y = w_vert[2].y;
  else if (w_vert[2].y > urb.y)
    urb.y = w_vert[2].y;

  if (w_vert[1].z < llf.z)
    llf.z = w_vert[1].z;
  else if (w_vert[1].z > urb.z)
    urb.z = w_vert[1].z;
  if (w_vert[2].z < llf.z)
    llf.z = w_vert[2].z;
  else if (w_vert[2].z > urb.z)
    urb.z = w_vert[2].z;
}


/***************************************************************************
distribute_wall:
  In: a wall 'w' that fully fits into partition 'p'
  Out: colliding_subparts - indices of all the subpartitions in a given partition
        where the wall is located
***************************************************************************/
// original name: distribute_wall
static void wall_subparts_collision_test(
    const Partition& p, const Wall& w,
    SubpartIndicesVector& colliding_subparts
) {
  vec3_t llf, urb; /* Bounding box for wall */
  float_t leeway = 1.0; /* Margin of error */

  vec3_t w_vert[VERTICES_IN_TRIANGLE];
  w_vert[0] = p.get_geometry_vertex(w.vertex_indices[0]);
  w_vert[1] = p.get_geometry_vertex(w.vertex_indices[1]);
  w_vert[2] = p.get_geometry_vertex(w.vertex_indices[2]);

  get_wall_bounding_box(w_vert, llf, urb);

  // min
  if (llf.x < -leeway)
    leeway = -llf.x;
  if (llf.y < -leeway)
    leeway = -llf.y;
  if (llf.z < -leeway)
    leeway = -llf.z;

  // max
  if (urb.x > leeway)
    leeway = urb.x;
  if (urb.y > leeway)
    leeway = urb.y;
  if (urb.z > leeway)
    leeway = urb.z;
  leeway = EPS + leeway * EPS;
  if (p.get_world_constants().use_expanded_list) {
    leeway += p.get_world_constants().rx_radius_3d;
  }

  llf = llf - vec3_t(leeway);
  urb = urb + vec3_t(leeway);

  // let's assume for now that we are placing a cube with corners llf and urb,
  // fing what are the min and max parition indices
  ivec3_t min_subpart_indices, max_subpart_indices;
  p.get_subpart_3d_indices(llf, min_subpart_indices);
  p.get_subpart_3d_indices(urb, max_subpart_indices);

  // do we fit into just one subparition?
  if (max_subpart_indices == min_subpart_indices) {
    colliding_subparts.push_back(
        p.get_subpartition_index_from_3d_indices(min_subpart_indices)
    );
    return;
  }

#if 0
  // don't know yet what is this doing, somehow it is trying to find some subvolume
  // but don't know which one

  // probably just an optimization

  // this is needed for the following code (if enabled, fix loop afterwards)
  max_subpart_indices += ivec3_t(1); // max is incremented by 1 in each dimension

  cent.x = 0.33333333333 * (w_vert[0].x + w_vert[1].x + w_vert[2].x);
  cent.y = 0.33333333333 * (w_vert[0].y + w_vert[1].y + w_vert[2].y);
  cent.z = 0.33333333333 * (w_vert[0].z + w_vert[1].z + w_vert[2].z);


  for (i = x_min; i < x_max; i++) {
    if (cent.x < world->x_partitions[i])
      break;
  }
  for (j = y_min; j < y_max; j++) {
    if (cent.y < world->y_partitions[j])
      break;
  }
  for (k = z_min; k < z_max; k++) {
    if (cent.z < world->z_partitions[k])
      break;
  }

  h = (k - 1) +
      (world->nz_parts - 1) * ((j - 1) + (world->ny_parts - 1) * (i - 1));
  where_am_i = localize_wall(w, world->subvol[h].local_storage);
  if (where_am_i == NULL)
    return NULL;
#endif

  // go through all subparts in llf-urb box and select only those that cross our wall
  for (int x = min_subpart_indices.x; x <= max_subpart_indices.x; x++) {
    for (int y = min_subpart_indices.y; y <= max_subpart_indices.y; y++) {
      for (int z = min_subpart_indices.z; z <= max_subpart_indices.z; z++) {

        vec3_t subpart_llf, subpart_urb;
        subpart_index_t subpart_index = p.get_subpart_index_from_3d_indices(x, y, z);
        p.get_subpart_llf_point(subpart_index, subpart_llf);
        p.get_subpart_urb_point_from_llf(subpart_llf, subpart_urb);

        if (WallUtil::wall_in_box(p, w, subpart_llf, subpart_urb) != 0) {
          colliding_subparts.push_back(subpart_index);
        }
      }
    }
  }
}



/***************************************************************************
find_edge_point:
  In: here: a wall
      loc: a point in the coordinate system of that wall where we are now
           (assumed to be on or inside triangle)
      disp: a 2D displacement vector to move
      edgept: a place to store the coordinate on the edge, if we hit it
  Out: index of the edge we hit (0, 1, or 2), or -1 if the new location
       is within the wall, or -2 if we can't tell.  If the result is
       0, 1, or 2, edgept is set to the new location.
***************************************************************************/
static int find_edge_point(
    const Wall& here,
    const vec2_t& loc,
    const vec2_t& disp,
    vec2_t& edgept
) {
  float_t lxd = determinant2(loc, disp);

  float_t lxc1 = -loc.v * here.uv_vert1_u;
  float_t dxc1 = -disp.v * here.uv_vert1_u;

  // Make sure that the displacement vector isn't on v0v1
  float_t f, s, t;
  if (dxc1 < -EPS || dxc1 > EPS) {
    f = 1.0 / dxc1; /* f>0 is passing outwards */
    s = -lxd * f;
    if (0.0 < s && s < 1.0 && f > 0.0) {
      t = -lxc1 * f;
      if (EPS < t && t < 1.0) {
        edgept = loc + vec2_t(t) * disp;
        return 0;
      }
      else if (t > 1.0 + EPS) {
        return -1;
      }
      /* else can't tell if we hit this edge, assume not */
    }
  }

  float_t  lxc2 = determinant2(loc, here.uv_vert2);
  float_t  dxc2 = determinant2(disp, here.uv_vert2);

  // Make sure that the displacement vector isn't on v1v2
  if (dxc2 < -EPS || dxc2 > EPS) {
    f = 1.0 / dxc2; /* f<0 is passing outwards */
    s = 1.0 + lxd * f;
    if (0.0 < s && s < 1.0 && f < 0.0) {
      t = -lxc2 * f;
      if (EPS < t && t < 1.0) {
        edgept = loc + vec2_t(t) * disp;
        return 2;
      }
      else if (t > 1.0 + EPS) {
        return -1;
      }
      /* else can't tell */
    }
  }

  f = dxc2 - dxc1;

  if (f < -EPS || f > EPS) {
    f = 1.0 / f; /* f>0 is passing outwards */
    s = -(lxd + dxc1) * f;
    if (0.0 < s && s < 1.0 && f > 0.0) {
      t = (here.uv_vert1_u * here.uv_vert2.v + lxc1 - lxc2) * f;
      if (EPS < t && t < 1.0) {
        edgept = loc + vec2_t(t) * disp;
        return 1;
      }
      else if (t > 1.0 + EPS) {
        return -1;
      }
      /* else can't tell */
    }
  }

  return -2; /* Couldn't tell whether we hit or not--calling function should
                pick another displacement */
}



/***************************************************************************
traverse_surface:
  In: here: a wall
      loc: a point in the coordinate system of that wall
      which: which edge to travel off of
      newloc: a vector to set for the new wall
  Out: NULL if the edge is not shared, or a pointer to the wall in that
       direction if it is shared. newloc is set to loc in the coordinate system
       of the new wall (after flattening the walls along their shared edge)
***************************************************************************/
static wall_index_t traverse_surface(const Wall& here, const vec2_t& loc, int which_edge,
  vec2_t& newloc) {

  wall_index_t there;

  const Edge& e = here.edges[which_edge];

  if (e.forward_index == here.index) {
    /* Apply forward transform to loc */
    there = e.backward_index;

    /* rotation */
    vec2_t tmp;
    tmp.u = e.cos_theta * loc.u + e.sin_theta * loc.v;
    tmp.v = -e.sin_theta * loc.u + e.cos_theta * loc.v;

    /* translation */
    newloc = tmp + e.translate;

    return there;
  }
  else {
    /* Apply inverse transform to loc */
    there = e.forward_index;

    /* inverse translation */
    vec2_t tmp;
    tmp = loc - e.translate;

    /* inverse rotation */
    newloc.u = e.cos_theta * tmp.u - e.sin_theta * tmp.v;
    newloc.v = e.sin_theta * tmp.u + e.cos_theta * tmp.v;

    return there;
  }

  return WALL_INDEX_INVALID;
}


/**************************************************************************
same_side:
        In: two points p1 and p2
            line defined by the points a and b
        Out: returns 1 if points p1 and p2 are on the same side of the line
             defined by the points a and b
**************************************************************************/
static int same_side(const vec3_t& p1, const vec3_t& p2, const vec3_t& a, const vec3_t& b) {
  vec3_t cp1, cp2, b_a, p1_a, p2_a;
  b_a = b - a;
  p1_a = p1 - a;
  p2_a = p2 - a;

  cp1 = cross(b_a, p1_a);
  cp2 = cross(b_a, p2_a);

  if (dot(cp1, cp2) >= 0) {
    return 1;
  }
  else {
    return 0;
  }
}


/************************************************************************
point_in_triangle:
        In: point p
            triangle defined by points a,b,c
        Out: returns 1 if point p is in the triangle defined by
             points (a,b,c) or lies on edges (a,b), (b,c) or (a,c).
        Note: If point p coincides with vertices (a,b,c) we consider that p
              is in the triangle.
************************************************************************/
static bool point_in_triangle(const vec3_t& p, const vec3_t& a, const vec3_t& b, const vec3_t& c) {

  if (same_side(p, a, b, c) && same_side(p, b, a, c) && same_side(p, c, a, b)) {
    return 1;
  }

  if (((!distinguishable_f(p.x, a.x, EPS)) &&
       (!distinguishable_f(p.y, a.y, EPS)) &&
       (!distinguishable_f(p.z, a.z, EPS))) ||
      ((!distinguishable_f(p.x, b.x, EPS)) &&
       (!distinguishable_f(p.y, b.y, EPS)) &&
       (!distinguishable_f(p.z, b.z, EPS))) ||
      ((!distinguishable_f(p.x, c.x, EPS)) &&
       (!distinguishable_f(p.y, c.y, EPS)) &&
       (!distinguishable_f(p.z, c.z, EPS)))) {
    return 1;
  }

  return 0;
}


/*******************************************************************
cross2D:
   In: 2D vectors a and b
   Out: 2D pseudo cross product Dot(Perp(a0,b)
   Note: The code adapted from "Real-Time Collision Detection" by
              Christer Ericson, p.205

*******************************************************************/
static float_t cross2D(const vec2_t& a, const vec2_t& b) {
  return ((a.v) * (b.u) - (a.u) * (b.v));
}

/*********************************************************************
point_in_triangle_2D:
   In: point p
       triangle defined by vertices a, b, c
   Out: Returns 1 if point p is inside the above defined triangle,
        and 0 otherwise.
        Note: The code adapted from "Real-Time Collision Detection" by
              Christer Ericson, p.206
***********************************************************************/
static bool point_in_triangle_2D(const vec2_t& p, const vec2_t& a,
                         const vec2_t& b, const vec2_t& c) {
  float_t pab, pbc, pca;

  pab = cross2D(p - a, b - a);
  pbc = cross2D(p - b, c - b);
  /* if P left of one of AB and BC and right of the other, not inside triangle
     - (pab and pbc have different signs */
  if (((pab > 0) && (pbc < 0)) || ((pab < 0) && (pbc > 0))) {
    return false;
  }

  pca = cross2D(p - c, a - c);
  /* if P left of one of AB and CA and right of the other, not inside triangle
   * - pab and pca have different signs */
  if (((pab > 0) && (pca < 0)) || ((pab < 0) && (pca > 0))) {
    return false;
  }

  /* if P left or right of all edges, so must be in (or on) the triangle */
  return true;
}


/***************************************************************************
closest_pt_point_triangle:
  In:  p - point
       a,b,c - vectors defining the vertices of the triangle.
  Out: final_result - closest point on triangle ABC to a point p.
       The code is adapted from "Real-time Collision Detection" by Christer
Ericson, ISBN 1-55860-732-3, p.141.

***************************************************************************/
static void closest_pt_point_triangle(
    const vec3_t& p, const vec3_t& a,
    const vec3_t& b, const vec3_t& c,
    vec3_t& final_result
) {

  /* Check if P in vertex region outside A */
  vec3_t ab, ac, ap;
  ab = b - a;
  ac = c - a;
  ap = p - a;

  float_t d1, d2;
  d1 = dot(ab, ap);
  d2 = dot(ac, ap);

  if (d1 <= 0.0 && d2 <= 0.0) {
    final_result = a;
    return;
  }

  /* Check if P in vertex region outside B */
  vec3_t bp = p - b;
  float_t d3, d4;
  d3 = dot(ab, bp);
  d4 = dot(ac, bp);
  if (d3 >= 0.0 && d4 <= d3) {
    final_result = b;
    return;
  }

  /* Check if P in edge region of AB, if so return projection of P onto AB */
  float_t v, w;
  vec3_t result1;
  float_t vc = d1 * d4 - d3 * d2;
  if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0) {
    v = d1 / (d1 - d3);
    result1 = ab * vec3_t(v);
    final_result = a + result1;
    return; /* barycentric coordinates (1-v,v,0) */
  }

  /* Check if P in vertex region outside C */
  vec3_t cp = p - c;
  float_t d5, d6;
  d5 = dot(ab, cp);
  d6 = dot(ac, cp);
  if (d6 >= 0.0 && d5 <= d6) {
    final_result = c;
    return;
  }

  /* Check if P in edge region of AC, if so return projection of P onto AC */
  float_t vb = d5 * d2 - d1 * d6;
  if (vb <= 0.0f && d2 >= 0.0f && d6 <= 0.0f) {
    w = d2 / (d2 - d6);
    result1 = ac * vec3_t(w);
    final_result = a + result1;
    return; /* barycentric coordinates (0, 1-w,w) */
  }

  /* Check if P in edge region of BC, if so return projection of P onto BC */
  float_t va = d3 * d6 - d5 * d4;
  if (va <= 0.0f && (d4 - d3) >= 0.0f && (d5 - d6) >= 0.0f) {
    w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
    result1 = (c - b) * vec3_t(w);
    final_result = b + result1;
    return; /*barycentric coordinates (0,1-w, w) */
  }

  /* P inside face region. Compute Q through its barycentric
     coordinates (u,v,w) */
  float_t denom = 1.0f / (va + vb + vc);
  v = vb * denom;
  w = vc * denom;
  ab = ab * vec3_t(v);
  ac = ac * vec3_t(w);
  result1 = ab + ac;
  final_result = a + result1;
  return; /* = u*a + v*b + w*c, u = va * denom = 1.0f - v -w */
}


/***************************************************************************
closest_interior_point:
  In: a point in 3D
      a wall
      the surface coordinates of the closest interior point on the wall
      how far away the point can be before we give up
  Out: return the distance^2 between the input point and closest point.
       Sets closest interior point.
  Note: the search distance currently isn't used.  This function is just
        a wrapper for closest_pt_point_triangle.  If the closest point is
        on an edge or corner, we scoot the point towards the centroid of
        the triangle so we're contained fully within the triangle.
***************************************************************************/

static float_t closest_interior_point(
    const Partition& p,
    const vec3_t& pt,
    const Wall& w,
    vec2_t& ip) {
#ifdef DEBUG_CLOSEST_INTERIOR_POINT
  std::cout << "closest_interior_point: " << pt << "\n";;
  w.dump(p, "", true);
#endif


  vec3_t v;
  const vec3_t& w_vert0 = p.get_wall_vertex(w, 0);
  const vec3_t& w_vert1 = p.get_wall_vertex(w, 1);
  const vec3_t& w_vert2 = p.get_wall_vertex(w, 2);

  closest_pt_point_triangle(pt, w_vert0, w_vert1, w_vert2, v);

  ip = xyz2uv(p, v, w);

  /* Check to see if we're lying on an edge; if so, scoot towards centroid. */
  /* ip lies on edge of wall if cross products are zero */

  int give_up_ctr = 0;
  int give_up = 10;
  float_t a1 = ip.u * w.uv_vert2.v - ip.v * w.uv_vert2.u;
  float_t a2 = w.uv_vert1_u * ip.v;
  vec2_t vert_0(0.0);
  vec2_t vert_1(w.uv_vert1_u, 0.0);

  while (
      give_up_ctr < give_up
      && (!distinguishable_f(ip.v, 0, EPS)
          || !distinguishable_f(a1, 0, EPS)
          || !distinguishable_f(a1 + a2, 2.0 * w.area, EPS)
          || !point_in_triangle_2D(ip, vert_0, vert_1, w.uv_vert2))
   ) {
    /* Move toward centroid. It's possible for this movement to be so small
     * that we are essentially stuck in this loop, so bail out after a set
     * number of tries. The number chosen is somewhat arbitrary. In most cases,
     * one try is sufficent. */
    ip.u = (1.0 - 5 * EPS) * ip.u +
            5 * EPS * 0.333333333333333 * (w.uv_vert1_u + w.uv_vert2.u);
    ip.v = (1.0 - 5 * EPS) * ip.v +
            5 * EPS * 0.333333333333333 * w.uv_vert2.v;

    a1 = ip.u * w.uv_vert2.v - ip.v * w.uv_vert2.u;
    a2 = w.uv_vert1_u * ip.v;

    give_up_ctr++;
  }

  float_t res = len3_squared(v - pt);

#ifdef DEBUG_CLOSEST_INTERIOR_POINT
  std::cout << "res: " << res << ", ip: " << ip << "\n";
#endif

  return res;
}


// returns the closest distance
// wall must belong to the same object and region
static float_t find_closest_wall(
    const Partition& p,
    const vec3_t& pos, const wall_index_t wall_that_moved_molecule,
    wall_index_t& best_wall_index,
    vec2_t& best_wall_pos2d
) {

  const Wall& moved_wall = p.get_wall(wall_that_moved_molecule);

  best_wall_index = WALL_INDEX_INVALID;
  float_t best_d2 = GIGANTIC4;
  best_wall_pos2d = vec2_t(0);

  // NOTE: use wall_that_moved_molecule to optimize the search,
  // however for now we are using the same approach as in mcell3 because
  // we need to match the behavior
  const GeometryObject& obj = p.get_geometry_object(moved_wall.object_index);

  for (const wall_index_t wall_index: obj.wall_indices) {
    const Wall& w = p.get_wall(wall_index);

    // TODO: regions - we are checking only object id here
    if (!moved_wall.is_same_region(w)) {
      continue;
    }

    vec2_t wall_pos2d;
    float_t d2 = closest_interior_point(p, pos, w, wall_pos2d);

    if (d2 <= best_d2) { // the <= is to emulate behavior of mcell3 that goes through the walls in opposite order
      best_d2 = d2;
      best_wall_index = w.index;
      best_wall_pos2d = wall_pos2d;
    }
  }

  return best_d2;
}


static bool is_point_above_plane_defined_by_wall(const Partition& p, const Wall& w, const vec3_t& pos) {

  // make vector pointing from any point to our position
  vec3_t w0_pos = pos - p.get_wall_vertex(w, 0);

  // dot product with normal gives ||a|| * ||b|| * cos(phi)
  float_t dot_prod = dot(w0_pos, w.normal);
  assert(!cmp_eq(dot_prod, 0, EPS) && "Checked point is on the plane");
  return dot_prod > 0;
}


} // namespace GeometryUtil

} // namespace MCell

#endif // SRC4_GEOMETRY_UTILS_INC_
