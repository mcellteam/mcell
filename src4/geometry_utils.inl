/******************************************************************************
 *
 * Copyright (C) 2006-2017 by
 * The Salk Institute for Biological Studies and
 * Pittsburgh Supercomputing Center, Carnegie Mellon University
 *
 * Use of this source code is governed by an MIT-style
 * license that can be found in the LICENSE file or at
 * https://opensource.org/licenses/MIT.
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

#include "wall_utils.inl"

namespace MCell {

namespace GeometryUtils {

static Vec2 xyz2uv(const Partition& p, const Vec3& a, const Wall& w) {
  Vec2 res;

  if (w.has_initialized_grid()) {
    const Grid& g = w.grid;

    res.u = a.x * w.unit_u.x + a.y * w.unit_u.y + a.z * w.unit_u.z -
           g.vert0.u;
    res.v = a.x * w.unit_v.x + a.y * w.unit_v.y + a.z * w.unit_v.z -
           g.vert0.v;
  }
  else {
    Vec3 pos = a - p.get_wall_vertex(w, 0);
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
    const Vec3 w_vert[VERTICES_IN_TRIANGLE],
    Vec3& llf, Vec3& urb
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
  Vec3 llf, urb; /* Bounding box for wall */
  pos_t leeway = 1; /* Margin of error */

  Vec3 w_vert[VERTICES_IN_TRIANGLE];
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
  leeway = POS_EPS + leeway * POS_EPS;
  if (p.config.use_expanded_list) {
    leeway += p.config.rxn_radius_3d;
  }

  Vec3 leeway3 = Vec3(leeway);
  llf = llf - leeway3;
  urb = urb + leeway3;

  // let's assume for now that we are placing a cube with corners llf and urb,
  // fing what are the min and max parition indices
  IVec3 min_subpart_indices, max_subpart_indices;
  p.get_subpart_3d_indices(llf, min_subpart_indices);
  p.get_subpart_3d_indices(urb, max_subpart_indices);

#if 0
  // don't know yet what is this doing, somehow it is trying to find some subvolume
  // but don't know which one

  // probably just an optimization

  // this is needed for the following code (if enabled, fix loop afterwards)
  max_subpart_indices += IVec3(1); // max is incremented by 1 in each dimension

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

        Vec3 subpart_llf, subpart_urb;
        subpart_index_t subpart_index = p.get_subpart_index_from_3d_indices(x, y, z);
        p.get_subpart_llf_point(subpart_index, subpart_llf);
        p.get_subpart_urb_point_from_llf(subpart_llf, subpart_urb);

        subpart_llf = subpart_llf - leeway3;
        subpart_urb = subpart_urb + leeway3;

        if (WallUtils::wall_in_box(p, w, subpart_llf, subpart_urb) != 0) {
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
  Out: index of the edge we hit (0, 1, or 2), or 3 if the new location
       is within the wall, or 4 if we can't tell.  If the result is
       0, 1, or 2, edgept is set to the new location.
***************************************************************************/
static edge_index_t find_edge_point(
    const Wall& here,
    const Vec2& loc,
    const Vec2& disp,
    Vec2& edgept
) {
  pos_t lxd = determinant2(loc, disp);

  pos_t lxc1 = -loc.v * here.uv_vert1_u;
  pos_t dxc1 = -disp.v * here.uv_vert1_u;

  // Make sure that the displacement vector isn't on v0v1
  pos_t f, s, t;
  if (dxc1 < -POS_EPS || dxc1 > POS_EPS) {
    f = 1 / dxc1; /* f>0 is passing outwards */
    s = -lxd * f;
    if (0 < s && s < 1 && f > 0) {
      t = -lxc1 * f;
      if (POS_EPS < t && t < 1) {
        edgept = loc + Vec2(t) * disp;
        return EDGE_INDEX_0;
      }
      else if (t > 1 + POS_EPS) {
        return EDGE_INDEX_WITHIN_WALL;
      }
      /* else can't tell if we hit this edge, assume not */
    }
  }

  pos_t lxc2 = determinant2(loc, here.uv_vert2);
  pos_t dxc2 = determinant2(disp, here.uv_vert2);

  // Make sure that the displacement vector isn't on v1v2
  if (dxc2 < -POS_EPS || dxc2 > POS_EPS) {
    f = 1 / dxc2; /* f<0 is passing outwards */
    s = 1 + lxd * f;
    if (0 < s && s < 1 && f < 0) {
      t = -lxc2 * f;
      if (POS_EPS < t && t < 1) {
        edgept = loc + Vec2(t) * disp;
        return EDGE_INDEX_2;
      }
      else if (t > 1 + POS_EPS) {
        return EDGE_INDEX_WITHIN_WALL;
      }
      /* else can't tell */
    }
  }

  f = dxc2 - dxc1;

  if (f < -POS_EPS || f > POS_EPS) {
    f = 1 / f; /* f>0 is passing outwards */
    s = -(lxd + dxc1) * f;
    if (0 < s && s < 1 && f > 0) {
      t = (here.uv_vert1_u * here.uv_vert2.v + lxc1 - lxc2) * f;
      if (POS_EPS < t && t < 1) {
        edgept = loc + Vec2(t) * disp;
        return EDGE_INDEX_1;
      }
      else if (t > 1 + POS_EPS) {
        return EDGE_INDEX_WITHIN_WALL;
      }
      /* else can't tell */
    }
  }

  return EDGE_INDEX_CANNOT_TELL; /* Couldn't tell whether we hit or not--calling function should
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
static wall_index_t traverse_surface(const Wall& here, const Vec2& loc, edge_index_t which_edge,
  Vec2& newloc) {

  wall_index_t there;

  const Edge& e = here.edges[which_edge];

  if (e.forward_index == here.index) {
    /* Apply forward transform to loc */
    there = e.backward_index;

    /* rotation */
    Vec2 tmp;
    tmp.u = e.get_cos_theta() * loc.u + e.get_sin_theta() * loc.v;
    tmp.v = -e.get_sin_theta() * loc.u + e.get_cos_theta() * loc.v;

    /* translation */
    newloc = tmp + e.get_translate();

    return there;
  }
  else {
    /* Apply inverse transform to loc */
    there = e.forward_index;

    /* inverse translation */
    Vec2 tmp;
    tmp = loc - e.get_translate();

    /* inverse rotation */
    newloc.u = e.get_cos_theta() * tmp.u - e.get_sin_theta() * tmp.v;
    newloc.v = e.get_sin_theta() * tmp.u + e.get_cos_theta() * tmp.v;

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
static int same_side(const Vec3& p1, const Vec3& p2, const Vec3& a, const Vec3& b) {
  Vec3 cp1, cp2, b_a, p1_a, p2_a;
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
static bool point_in_triangle(const Vec3& p, const Vec3& a, const Vec3& b, const Vec3& c) {

  if (same_side(p, a, b, c) && same_side(p, b, a, c) && same_side(p, c, a, b)) {
    return 1;
  }

  if (((!distinguishable_p(p.x, a.x, POS_EPS)) &&
       (!distinguishable_p(p.y, a.y, POS_EPS)) &&
       (!distinguishable_p(p.z, a.z, POS_EPS))) ||
      ((!distinguishable_p(p.x, b.x, POS_EPS)) &&
       (!distinguishable_p(p.y, b.y, POS_EPS)) &&
       (!distinguishable_p(p.z, b.z, POS_EPS))) ||
      ((!distinguishable_p(p.x, c.x, POS_EPS)) &&
       (!distinguishable_p(p.y, c.y, POS_EPS)) &&
       (!distinguishable_p(p.z, c.z, POS_EPS)))) {
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
static pos_t cross2D(const Vec2& a, const Vec2& b) {
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
static bool point_in_triangle_2D(const Vec2& p, const Vec2& a,
                         const Vec2& b, const Vec2& c) {
  pos_t pab, pbc, pca;

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
    const Vec3& p, const Vec3& a,
    const Vec3& b, const Vec3& c,
    Vec3& final_result
) {

  /* Check if P in vertex region outside A */
  Vec3 ab, ac, ap;
  ab = b - a;
  ac = c - a;
  ap = p - a;

  pos_t d1, d2;
  d1 = dot(ab, ap);
  d2 = dot(ac, ap);

  if (d1 <= 0 && d2 <= 0) {
    final_result = a;
    return;
  }

  /* Check if P in vertex region outside B */
  Vec3 bp = p - b;
  pos_t d3, d4;
  d3 = dot(ab, bp);
  d4 = dot(ac, bp);
  if (d3 >= 0 && d4 <= d3) {
    final_result = b;
    return;
  }

  /* Check if P in edge region of AB, if so return projection of P onto AB */
  pos_t v, w;
  Vec3 result1;
  pos_t vc = d1 * d4 - d3 * d2;
  if (vc <= 0 && d1 >= 0 && d3 <= 0) {
    v = d1 / (d1 - d3);
    result1 = ab * Vec3(v);
    final_result = a + result1;
    return; /* barycentric coordinates (1-v,v,0) */
  }

  /* Check if P in vertex region outside C */
  Vec3 cp = p - c;
  pos_t d5, d6;
  d5 = dot(ab, cp);
  d6 = dot(ac, cp);
  if (d6 >= 0 && d5 <= d6) {
    final_result = c;
    return;
  }

  /* Check if P in edge region of AC, if so return projection of P onto AC */
  pos_t vb = d5 * d2 - d1 * d6;
  if (vb <= 0 && d2 >= 0 && d6 <= 0) {
    w = d2 / (d2 - d6);
    result1 = ac * Vec3(w);
    final_result = a + result1;
    return; /* barycentric coordinates (0, 1-w,w) */
  }

  /* Check if P in edge region of BC, if so return projection of P onto BC */
  pos_t va = d3 * d6 - d5 * d4;
  if (va <= 0 && (d4 - d3) >= 0 && (d5 - d6) >= 0) {
    w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
    result1 = (c - b) * Vec3(w);
    final_result = b + result1;
    return; /*barycentric coordinates (0,1-w, w) */
  }

  /* P inside face region. Compute Q through its barycentric
     coordinates (u,v,w) */
  pos_t denom = 1 / (va + vb + vc);
  v = vb * denom;
  w = vc * denom;
  ab = ab * Vec3(v);
  ac = ac * Vec3(w);
  result1 = ab + ac;
  final_result = a + result1;
  return; /* = u*a + v*b + w*c, u = va * denom = 1 - v -w */
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

static pos_t closest_interior_point(
    const Partition& p,
    const Vec3& pt,
    const Wall& w,
    Vec2& ip) {
#ifdef DEBUG_CLOSEST_INTERIOR_POINT
  std::cout << "closest_interior_point: " << pt << "\n";;
  w.dump(p, "", true);
#endif

  Vec3 v;
  const Vec3& w_vert0 = p.get_wall_vertex(w, 0);
  const Vec3& w_vert1 = p.get_wall_vertex(w, 1);
  const Vec3& w_vert2 = p.get_wall_vertex(w, 2);

  closest_pt_point_triangle(pt, w_vert0, w_vert1, w_vert2, v);

  ip = xyz2uv(p, v, w);

  /* Check to see if we're lying on an edge; if so, scoot towards centroid. */
  /* ip lies on edge of wall if cross products are zero */

  int give_up_ctr = 0;
  int give_up = 10;
  pos_t a1 = ip.u * w.uv_vert2.v - ip.v * w.uv_vert2.u;
  pos_t a2 = w.uv_vert1_u * ip.v;
  Vec2 vert_0(0);
  Vec2 vert_1(w.uv_vert1_u, 0);

  while (
      give_up_ctr < give_up
      && (!distinguishable_p(ip.v, 0, POS_EPS) ||
          !distinguishable_p(a1, 0, POS_EPS) ||
          !distinguishable_p(a1 + a2, 2 * w.area, POS_EPS) ||
          !point_in_triangle_2D(ip, vert_0, vert_1, w.uv_vert2))
   ) {
    /* Move toward centroid. It's possible for this movement to be so small
     * that we are essentially stuck in this loop, so bail out after a set
     * number of tries. The number chosen is somewhat arbitrary. In most cases,
     * one try is sufficent. */
    ip.u = (1 - 5 * POS_EPS) * ip.u +
            5 * POS_EPS * 0.333333333333333 * (w.uv_vert1_u + w.uv_vert2.u);
    ip.v = (1 - 5 * POS_EPS) * ip.v +
            5 * POS_EPS * 0.333333333333333 * w.uv_vert2.v;

    a1 = ip.u * w.uv_vert2.v - ip.v * w.uv_vert2.u;
    a2 = w.uv_vert1_u * ip.v;

    give_up_ctr++;
  }

  pos_t res = len3_squared(v - pt);

#ifdef DEBUG_CLOSEST_INTERIOR_POINT
  std::cout << "res: " << res << ", ip: " << ip << "\n";
#endif

  return res;
}


static bool is_point_above_plane_defined_by_wall(const Partition& p, const Wall& w, const Vec3& pos) {

  // make vector pointing from any point to our position
  Vec3 w0_pos = pos - p.get_wall_vertex(w, 0);

  // dot product with normal gives ||a|| * ||b|| * cos(phi)
  pos_t dot_prod = dot(w0_pos, w.normal);
  assert(!cmp_eq(dot_prod, (pos_t)0.0, POS_EPS) && "Checked point is on the plane");
  return dot_prod > 0;
}


} // namespace GeometryUtil

} // namespace MCell

#endif // SRC4_GEOMETRY_UTILS_INC_
