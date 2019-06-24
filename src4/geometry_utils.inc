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

#include "diffuse_react_event.h"
#include "defines.h"
#include "world.h"
#include "partition.h"
#include "geometry.h"
#include "debug_config.h"

namespace mcell {

namespace geom_util {

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
    const wall_t& here,
    const vec2_t& loc,
    const vec2_t& disp,
    vec2_t& edgept
) {
  // TODO: cleanup

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
static wall_index_t traverse_surface(wall_t& here, const vec2_t& loc, int which_edge,
  vec2_t& newloc) {

  wall_index_t there;

  const edge_t& e = here.edges[which_edge];

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



/*******************************************************************
cross2D:
   In: 2D vectors a and b
   Out: 2D pseudo cross product Dot(Perp(a0,b)
   Note: The code adapted from "Real-Time Collision Detection" by
              Christer Ericson, p.205

*******************************************************************/
// FIXME - replace with determinant2
static float_t cross2D(const vec2_t& a, const vec2_t& b) {
  return ((a.v) * (b.u) - (a.u) * (b.v));
}

/*************************************************************************
vectorize2D:
   In: 2D vectors p1 and p2
   Out: Subtracts vector p1 from p2 and places result into p3
*************************************************************************/
static void vectorize2D(const vec2_t& p1, const vec2_t& p2, vec2_t& p3) {
  p3.u = p2.u - p1.u;
  p3.v = p2.v - p1.v;
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
static int point_in_triangle_2D(const vec2_t& p, const vec2_t& a,
                         const vec2_t& b, const vec2_t& c) {
  vec2_t p_minus_a, b_minus_a, p_minus_b, c_minus_b, p_minus_c,
      a_minus_c;
  float_t pab, pbc, pca;

  vectorize2D(a, p, p_minus_a);
  vectorize2D(a, b, b_minus_a);
  vectorize2D(b, p, p_minus_b);
  vectorize2D(b, c, c_minus_b);
  vectorize2D(c, p, p_minus_c);
  vectorize2D(c, a, a_minus_c);

  pab = cross2D(p_minus_a, b_minus_a);
  pbc = cross2D(p_minus_b, c_minus_b);
  /* if P left of one of AB and BC and right of the other, not inside triangle
     - (pab and pbc have different signs */
  if (((pab > 0) && (pbc < 0)) || ((pab < 0) && (pbc > 0))) {
    return 0;
  }

  pca = cross2D(p_minus_c, a_minus_c);
  /* if P left of one of AB and CA and right of the other, not inside triangle
   * - pab and pca have different signs */
  if (((pab > 0) && (pca < 0)) || ((pab < 0) && (pca > 0))) {
    return 0;
  }

  /* if P left or right of all edges, so must be in (or on) the triangle */
  return 1;
}



static int uv2grid(
    const vec2_t& v,
    const wall_t& w,
    const grid_t& g
) {
  float_t i, j;
  float_t u0, u1_u0;
  float_t striploc, striprem, stripeloc, striperem;
  int strip, stripe, flip, idx;
  vec2_t vert_0, vert_1;
  int tile_idx_0, tile_idx_mid, tile_idx_last;

  if (g.num_tiles == 1) {
    return 0;
  }

  /* find tile indices of the corner tiles */
  tile_idx_0 = 0;
  /* see function "move_strip_up()" */
  tile_idx_mid = g.num_tiles - 2 * g.num_tiles_along_axis + 1;
  tile_idx_last = g.num_tiles - 1;

  vert_0.u = vert_0.v = 0;
  vert_1.u = w.uv_vert1_u;
  vert_1.v = 0;

  if (!distinguishable_vec2(v, vert_0, EPS)) {
    return tile_idx_mid;
  }
  if (!distinguishable_vec2(v, vert_1, EPS)) {
    return tile_idx_0;
  }
  if (!distinguishable_vec2(v, w.uv_vert2, EPS)) {
    return tile_idx_last;
  }

  // TODO: just for debug?
  if (!(point_in_triangle_2D(v, vert_0, vert_1, w.uv_vert2))) {
    mcell_internal_error("Error in function 'uv2grid()': point is outside wall.");
  }

  i = v.u;
  j = v.v;

  striploc = j * g.inv_strip_wid;
  strip = (int)striploc;
  striprem = striploc - strip;

  strip = g.num_tiles_along_axis - strip - 1;

  u0 = j * g.vert2_slope;
  u1_u0 = w.uv_vert1_u - j * g.fullslope;

  stripeloc = ((i - u0) / u1_u0) * (((double)strip) + (1.0 - striprem));
  stripe = (int)(stripeloc);
  striperem = stripeloc - stripe;

  flip = (striperem < 1.0 - striprem) ? 0 : 1;
  idx = strip * strip + 2 * stripe + flip;

  if ((u_int)idx >= g.num_tiles) {
    mcell_internal_error(
        "Error in function 'xyz2grid()': returning tile index "
        "%d while wall has %u tiles",
        idx, g.num_tiles
    );
  }

  return idx;
}

} // namespace geom_util

} // namespace mcell

#endif // SRC4_GEOMETRY_UTILS_INC_
