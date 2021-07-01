/******************************************************************************
 *
 * Copyright (C) 2006-2017,2020 by
 * The Salk Institute for Biological Studies and
 * Pittsburgh Supercomputing Center, Carnegie Mellon University
 *
 * Use of this source code is governed by an MIT-style
 * license that can be found in the LICENSE file or at
 * https://opensource.org/licenses/MIT.
 *
******************************************************************************/

#include "grid_position.h"
#include "geometry.h"
#include "logging.h"

#include "geometry_utils.inl"
#include "grid_utils.inl"

namespace MCell {
namespace GridPosition {

using GeometryUtils::xyz2uv;
using GeometryUtils::uv2xyz;

/*************************************************************************
get_tile_vertices:
   In: surface grid
       tile index
       tile flip information (return value)
       first tile vertex R   (return value)
       second tile vertex S  (return value)
       third tile vertex T   (return value)
   Out: the tile vertices (R,S,T) coordinates are defined
   Note: the vertices (R,S,T) are listed clockwise.  For the upright tile
         (orientation similar to the wall) two vertices R and S are
         on the same line PQ parallel and closest to the u-axis.
         For the inverted tile (orientation inverted compared to the wall)
         two vertices S and T are on the same line XY parallel and furthest
         to the u-axis.

*************************************************************************/
static void get_tile_vertices(
    const Partition& p, const Grid& sg, const tile_index_t idx,
    int& flp,
    Vec2& R, Vec2& S, Vec2& T) {

  /* indices of the barycentric subdivision */
  int strip, stripe, flip;
  int root, rootrem;
  Vec2 P, Q, X, Y;
  /* length of the segments PQ and XY (see below) */
  pos_t pq, xy;
  /* cotangent of the angle formed between the u-axis
     and the wall edge opposite to the origin of the
     uv-coordinate system */
  pos_t cot_angle;

  /* Calculate strip, stripe, and flip indices from idx */
  root = (int)sqrt_p(idx);
  rootrem = idx - root * root;
  strip = sg.num_tiles_along_axis - root - 1;
  stripe = rootrem / 2;
  flip = rootrem - 2 * stripe;

  /* Let PQ to be the segment on the grid containing the vertex R
     and the one closest to the u-axis.  Let XY to be the segment
     containing the vertex T and the one furthest to the u-axis.
     Let point P to be on the left side of point Q.
     Let point X to be on the left side of point Y.
  */

  /* find v-coordinates of P, Q, X, Y */
  P.v = Q.v = strip / sg.strip_width_rcp;
  X.v = Y.v = (strip + 1) / sg.strip_width_rcp;

  /* find u-coordinates of P, Q, X, Y */
  P.u = (P.v) * (sg.vert2_slope);
  X.u = (X.v) * (sg.vert2_slope);

  const Wall& surface = p.get_wall(sg.wall_index);

  cot_angle = (surface.uv_vert1_u - surface.uv_vert2.u) /
              (surface.uv_vert2.v);
  Q.u = surface.uv_vert1_u - (Q.v) * cot_angle;
  Y.u = surface.uv_vert1_u - (Y.v) * cot_angle;

  pq = Q.u - P.u;
  if (idx == 0) {
    xy = 0;
  } else {
    xy = Y.u - X.u;
  }

  /* find coordinates of the tile vertices */
  /* For the upright tile the vertices R and S lie on
     the line PQ and vertex T - on the line XY.
     For the inverted tile only the vertex R lies on the line PQ
     while the vertices S and T lie on the line XY */

  if (flip == 1) {
    /* inverted tile */
    R.v = P.v;
    /* note: pq/(sg.num_tiles_along_axis - strip) tells us
             the number of slots on the line PQ */
    R.u = P.u + pq * (stripe + 1) / (sg.num_tiles_along_axis - strip);

    S.v = T.v = X.v;
    T.u = X.u + xy * stripe / (sg.num_tiles_along_axis - strip - 1);
    S.u = X.u + xy * (stripe + 1) / (sg.num_tiles_along_axis - strip - 1);

  } else {
    /* upright tile */
    R.v = S.v = P.v;
    T.v = X.v;

    /* note: pq/(sg.num_tiles_along_axis - strip) tells us
             the number of slots on the line PQ */
    R.u = P.u + pq * stripe / (sg.num_tiles_along_axis - strip);
    S.u = P.u + pq * (stripe + 1) / (sg.num_tiles_along_axis - strip);
    if (idx == 0) {
      T.u = X.u;
    } else {
      T.u = X.u + xy * stripe / (sg.num_tiles_along_axis - strip - 1);
    }
  }

  /* set the tile flip value */
  flp = flip;
}


/****************************************************************************
find_wall_vertex_for_corner_tile:
   In: surface grid
       tile index on the grid
   Out: corner tile should share one of it's vertices with the wall.
        Returns the shared wall vertex id (index in the array wall_vert[])
****************************************************************************/
int find_wall_vertex_for_corner_tile(const Grid& grid, const tile_index_t idx) {
  /* index of the wall vertex that is shared with the tile vertex */
  int vertex_id = 0;

  if (!GridUtils::is_corner_tile(grid, idx))
    mcell_internal_error("Function 'find_wall_vertex_for_corner_tile()' is "
                         "called for the tile that is not the corner tile.");

  if ((u_int)idx == (grid.num_tiles - 2 * (grid.num_tiles_along_axis) + 1)) {
    vertex_id = 0;
  }
  else if ((u_int)idx == (grid.num_tiles - 1)) {
    vertex_id = 1;
  }
  else if (idx == 0) {
    vertex_id = 2;
  }
  else {
    mcell_internal_error("Function 'find_wall_vertex_for_corner_tile()' is "
                         "called for the tile that is not the corner tile.");
  }

  return vertex_id;
}


/*****************************************************************************
get_product_shared_segment_pos:
   In: segment defined by points R_shared and S_shared
       point T
       position of the product (return value)
       ratios k1 and k2
   Out: coordinates of the point "prod" dividing the distance between T
       and midpont on (R_shared, S_shared) in the ratio k1:k2, so that
       midpoint_prod/midpoint_T = k1/k2
   Note: all points are on the plane
******************************************************************************/
static Vec2 get_product_shared_segment_pos(
    const Vec2& R_shared, const Vec2& S_shared, const Vec2& T, pos_t k1, pos_t k2) {

  Vec2 M; /* midpoint */

  /* find midpoint on the segment (R_shared - S_shared) */
  M.u = 0.5 * (R_shared.u + S_shared.u);
  M.v = 0.5 * (R_shared.v + S_shared.v);

  /* find coordinates of the  point such that internally
      divides the segment TM in the ratio (k1:k2)
      We want to place the product close to the common edge.
  */
  return Vec2(
      (k1 * T.u + k2 * M.u) / (k1 + k2),
      (k1 * T.v + k2 * M.v) / (k1 + k2));
}


/*****************************************************************************
get_product_shared_vertex_pos:
   In: point R_shared
       segment defined by points S and T
       position of the product (return value)
       ratios k1 and k2
   Out: coordinates of the point "prod" dividing the distance between R_shared
       and midpont on (S, T) in the ratio k1:k2, so that
       R_shared_prod/prod_midpoint = k1/k2
   Note: all points are on the plane
******************************************************************************/
static Vec2 get_product_shared_vertex_pos(
    const Vec2& R_shared, const Vec2& S, const Vec2& T, const pos_t k1, const pos_t k2) {

  struct vector2 M; /* midpoint */

  /* find midpoint on the segment ST */
  M.u = (pos_t)0.5 * (S.u + T.u);
  M.v = (pos_t)0.5 * (S.v + T.v);

  /* find coordinates of the  point such that internally
     divides the segment RM in the ratio (k1:k2)
  */
  return Vec2(
      (k1 * M.u + k2 * R_shared.u) / (k1 + k2),
      (k1 * M.v + k2 * R_shared.v) / (k1 + k2));
}


/*****************************************************************************
get_product_close_to_segment_endpoint_pos:
   In: segment defined by points S (start) and E (end)
       position of the product (return value)
       ratios k1 and k2
   Out: coordinates of the point "prod" dividing the distance between S and E
        in the ratio k1:k2, so that E_prod/S_prod = k1/k2
   Note: all points are on the plane
******************************************************************************/
static Vec2 get_product_close_to_segment_endpoint_pos(
    const Vec2& S, const Vec2& E, const pos_t k1, const pos_t k2) {

  return Vec2(
      (k1 * S.u + k2 * E.u) / (k1 + k2),
      (k1 * S.v + k2 * E.v) / (k1 + k2));
}


/*****************************************************************************
intersect_point_segment:
   In: point P and segment AB
   Out: 1 if the point P lies on the segment AB, and 0 - otherwise
******************************************************************************/
static bool intersect_point_segment(
    const Vec3& P, const Vec3& A, const Vec3& B) {

  Vec3 ba, pa;
  pos_t t;                    /* parameter in the line parametrical equation */
  pos_t ba_length, pa_length; /* length of the vectors */
  pos_t cosine_angle;         /* cosine of the angle between ba and pa */

  /* check for the end points */
  if (!distinguishable_vec3(P, A, POS_EPS)) {
    return true;
  }
  if (!distinguishable_vec3(P, B, POS_EPS)) {
    return true;
  }

  ba = B - A;
  pa = P - A;

  ba_length = len3(ba);
  pa_length = len3(pa);

  /* if point intersects segment, vectors pa and ba should be collinear */
  cosine_angle = dot(ba, pa) / (ba_length * pa_length);
  if (distinguishable_p(cosine_angle, 1, POS_EPS)) {
    return false;
  }

  /* Project P on AB, computing parameterized position d(t) = A + t(B - A) */
  t = dot(pa, ba) / dot(ba, ba);

  if ((t > 0) && (t < 1)) {
    return true;
  }

  return false;
}


/****************************************************************************
parallel_segments:
   In: segment defined by endpoints A, B
       segment defined by endpoints R, S
   Out: 1, if the segments are parallel.
        0, otherwise
*****************************************************************************/
bool parallel_segments(
    const Vec3& A, const Vec3& B,
    const Vec3& R, const Vec3& S) {

  double length;
  Vec3 prod; /* cross product */
  Vec3 ba, rs;

  ba = B - A;
  rs = R - S;
  prod = cross(ba, rs);
  length = len3(prod);

  if (!distinguishable_p(length, 0, POS_EPS)) {
    return true;
  }
  else {
    return false;
  }
}


/************************************************************************
find_closest_position:
  In: surface grid of the first tile
      first tile index
      surface grid of the second tile
      second (neighbor) tile index
  Out: position of the product on the first tile that is closest
       to the second tile. If the neighbor tiles have common edge
       this position happens to be very close to the center of the
       common edge but inward the first tile.  If the neighbor tiles
       have common vertex, this position happens to be very close to
       to the vertex but inward the first tile.
*************************************************************************/
// TODO: this function needs refactor
Vec2 find_closest_position(const Partition& p, const GridPos& gp1, const GridPos& gp2) {

  assert(gp1.is_assigned());
  assert(gp2.is_assigned());
  assert(!gp1.has_same_wall_and_grid(gp2));

  /* vertices of the first tile */
  Vec2 R, S, T;
  Vec3 R_3d, S_3d, T_3d;

  /* vertices of the second tile */
  Vec2 A, B, C;
  Vec3 A_3d, B_3d, C_3d;
  /* vertices A,B,C in the coordinate system RST */
  Vec2 A_new, B_new, C_new;

  /* the ratios in which we divide the segment */
  pos_t k1 = 1e-10; /* this is our good faith assumption */
  pos_t k2 = 1;

  int flip1; /* flip information about first tile */
  int flip2; /* flip information about second tile */

  int num_exact_shared_vertices = 0;
  /* flags */
  int R_shared = 0, S_shared = 0, T_shared = 0;
  int A_shared = 0, B_shared = 0, C_shared = 0;

  /* find out the vertices of the first tile where we will put the product */
  const Wall& wall1 = p.get_wall(gp1.wall_index);
  const Vec3& wall1_vert0 = p.get_wall_vertex(wall1, 0);
  const Grid& grid1 = wall1.grid;
  tile_index_t idx1 = gp1.tile_index;

  get_tile_vertices(p, grid1, idx1, flip1, R, S, T);

  /* the code below tries to increase accuracy for the corner tiles */
  if (GridUtils::is_corner_tile(grid1, idx1)) {
    /* find out the shared vertex */
    int shared_wall_vertex_id_1 = find_wall_vertex_for_corner_tile(grid1, idx1);
    /* note that vertices R, S, T followed clockwise rule */
    if (idx1 == 0) {
      T_3d = p.get_wall_vertex(wall1, shared_wall_vertex_id_1);
      R_3d = uv2xyz(R, wall1, wall1_vert0);
      S_3d = uv2xyz(S, wall1, wall1_vert0);
    }
    else if (idx1 == (grid1.num_tiles - 2 * (grid1.num_tiles_along_axis) + 1)) {
      R_3d = p.get_wall_vertex(wall1, shared_wall_vertex_id_1);
      S_3d = uv2xyz(S, wall1, wall1_vert0);
      T_3d = uv2xyz(T, wall1, wall1_vert0);
    }
    else {
      S_3d = p.get_wall_vertex(wall1, shared_wall_vertex_id_1);
      R_3d = uv2xyz(R, wall1, wall1_vert0);
      T_3d = uv2xyz(T, wall1, wall1_vert0);
    }
  }
  else {
    R_3d = uv2xyz(R, wall1, wall1_vert0);
    S_3d = uv2xyz(S, wall1, wall1_vert0);
    T_3d = uv2xyz(T, wall1, wall1_vert0);
  }

  const Wall& wall2 = p.get_wall(gp2.wall_index);
  const Vec3& wall2_vert0 = p.get_wall_vertex(wall2, 0);
  const Grid& grid2 = wall2.grid;
  tile_index_t idx2 = gp2.tile_index;

  get_tile_vertices(p, grid2, idx2, flip1, A, B, C);

  /* the code below tries to increase accuracy for the corner tiles */
  // TODO: make a function out of this, same code as above
  if (GridUtils::is_corner_tile(grid2, idx2)) {
    /* find out the shared vertex */
    int shared_wall_vertex_id_1 = find_wall_vertex_for_corner_tile(grid2, idx2);
    /* note that vertices A, B, C followed clockwise rule */
    if (idx2 == 0) {
      C_3d = p.get_wall_vertex(wall2, shared_wall_vertex_id_1);
      A_3d = uv2xyz(A, wall2, wall2_vert0);
      B_3d = uv2xyz(B, wall2, wall2_vert0);
    }
    else if (idx2 == (grid2.num_tiles - 2 * (grid2.num_tiles_along_axis) + 1)) {
      A_3d = p.get_wall_vertex(wall2, shared_wall_vertex_id_1);
      B_3d = uv2xyz(B, wall2, wall2_vert0);
      C_3d = uv2xyz(C, wall2, wall2_vert0);
    }
    else {
      B_3d = p.get_wall_vertex(wall2, shared_wall_vertex_id_1);
      A_3d = uv2xyz(A, wall2, wall2_vert0);
      C_3d = uv2xyz(C, wall2, wall2_vert0);
    }
  }
  else {
    A_3d = uv2xyz(A, wall2, wall2_vert0);
    B_3d = uv2xyz(B, wall2, wall2_vert0);
    C_3d = uv2xyz(C, wall2, wall2_vert0);
  }

  /* find shared vertices */
  if (gp1.wall_index == gp2.wall_index) {
    if (!distinguishable_vec2(R, A, POS_EPS) ||
        (!distinguishable_vec2(R, B, POS_EPS)) ||
        (!distinguishable_vec2(R, C, POS_EPS))) {
      num_exact_shared_vertices++;
      R_shared = 1;
    }
    if (!distinguishable_vec2(S, A, POS_EPS) ||
        (!distinguishable_vec2(S, B, POS_EPS)) ||
        (!distinguishable_vec2(S, C, POS_EPS))) {
      num_exact_shared_vertices++;
      S_shared = 1;
    }
    if (!distinguishable_vec2(T, A, POS_EPS) ||
        (!distinguishable_vec2(T, B, POS_EPS)) ||
        (!distinguishable_vec2(T, C, POS_EPS))) {
      num_exact_shared_vertices++;
      T_shared = 1;
    }

  } else {
    /* below there are cases when the grid structures on the neighbor
       walls are not shifted relative to one another */
    if (!distinguishable_vec3(R_3d, A_3d, POS_EPS) ||
        (!distinguishable_vec3(R_3d, B_3d, POS_EPS)) ||
        (!distinguishable_vec3(R_3d, C_3d, POS_EPS))) {
      num_exact_shared_vertices++;
      R_shared = 1;
    }

    if (!distinguishable_vec3(S_3d, A_3d, POS_EPS) ||
        (!distinguishable_vec3(S_3d, B_3d, POS_EPS)) ||
        (!distinguishable_vec3(S_3d, C_3d, POS_EPS))) {
      num_exact_shared_vertices++;
      S_shared = 1;
    }

    if (!distinguishable_vec3(T_3d, A_3d, POS_EPS) ||
        (!distinguishable_vec3(T_3d, B_3d, POS_EPS)) ||
        (!distinguishable_vec3(T_3d, C_3d, POS_EPS))) {
      num_exact_shared_vertices++;
      T_shared = 1;
    }
  }

  if (num_exact_shared_vertices == 1) {
    if (R_shared) {
      return get_product_shared_vertex_pos(R, S, T, k1, k2);
    }
    else if (S_shared) {
      return get_product_shared_vertex_pos(S, R, T, k1, k2);
    }
    else { /*T is shared */
      return get_product_shared_vertex_pos(T, R, S, k1, k2);
    }
  }

  if (num_exact_shared_vertices == 2) {
    if (R_shared && S_shared) {
      return get_product_shared_segment_pos(R, S, T, k1, k2);
    }
    else if (R_shared && T_shared) {
      return get_product_shared_segment_pos(R, T, S, k1, k2);
    }
    else { /*S_shared and T_shared */
      return get_product_shared_segment_pos(S, T, R, k1, k2);
    }
  }

  if (num_exact_shared_vertices == 0) {
    /* below are the cases when the grids on the neighbor walls
       are shifted relative to one another */
    /* find out whether the vertices of one tile cross the sides of
       another tile */

    if ((intersect_point_segment(S_3d, A_3d, B_3d)) ||
        (intersect_point_segment(S_3d, B_3d, C_3d)) ||
        (intersect_point_segment(S_3d, A_3d, C_3d))) {
      S_shared = 1;
    }

    if ((intersect_point_segment(R_3d, A_3d, B_3d)) ||
        (intersect_point_segment(R_3d, B_3d, C_3d)) ||
        (intersect_point_segment(R_3d, A_3d, C_3d))) {
      R_shared = 1;
    }

    if ((intersect_point_segment(T_3d, A_3d, B_3d)) ||
        (intersect_point_segment(T_3d, B_3d, C_3d)) ||
        (intersect_point_segment(T_3d, A_3d, C_3d))) {
      T_shared = 1;
    }

    if ((intersect_point_segment(A_3d, R_3d, S_3d)) ||
        (intersect_point_segment(A_3d, S_3d, T_3d)) ||
        (intersect_point_segment(A_3d, R_3d, T_3d))) {
      A_shared = 1;
    }

    if ((intersect_point_segment(B_3d, R_3d, S_3d)) ||
        (intersect_point_segment(B_3d, S_3d, T_3d)) ||
        (intersect_point_segment(B_3d, R_3d, T_3d))) {
      B_shared = 1;
    }

    if ((intersect_point_segment(C_3d, R_3d, S_3d)) ||
        (intersect_point_segment(C_3d, S_3d, T_3d)) ||
        (intersect_point_segment(C_3d, R_3d, T_3d))) {
      C_shared = 1;
    }

    /* two vertices shared from the same tile */
    if (R_shared && S_shared) {
      return get_product_shared_segment_pos(R, S, T, k1, k2);
    }
    else if (R_shared && T_shared) {
      return get_product_shared_segment_pos(R, T, S, k1, k2);
    }
    else if (S_shared && T_shared) {
      return get_product_shared_segment_pos(S, T, R, k1, k2);
    }

    /* two vertices shared from the same tile */
    if (A_shared && B_shared) {
      if (parallel_segments(A_3d, B_3d, R_3d, S_3d)) {
        A_new = xyz2uv(p, A_3d, wall1);
        B_new = xyz2uv(p, B_3d, wall1);
        return get_product_shared_segment_pos(A_new, B_new, T, k1, k2);
      }
      else if (parallel_segments(A_3d, B_3d, R_3d, T_3d)) {
        A_new = xyz2uv(p, A_3d, wall1);
        B_new = xyz2uv(p, B_3d, wall1);
        return get_product_shared_segment_pos(A_new, B_new, S, k1, k2);
      }
      else if (parallel_segments(A_3d, B_3d, S_3d, T_3d)) {
        A_new = xyz2uv(p, A_3d, wall1);
        B_new = xyz2uv(p, B_3d, wall1);
        return get_product_shared_segment_pos(A_new, B_new, R, k1, k2);
      }

    }
    else if (A_shared && C_shared) {
      if (parallel_segments(A_3d, C_3d, R_3d, S_3d)) {
        A_new = xyz2uv(p, A_3d, wall1);
        C_new = xyz2uv(p, C_3d, wall1);
        return get_product_shared_segment_pos(A_new, C_new, T, k1, k2);
      }
      else if (parallel_segments(A_3d, C_3d, R_3d, T_3d)) {
        A_new = xyz2uv(p, A_3d, wall1);
        C_new = xyz2uv(p, C_3d, wall1);
        return get_product_shared_segment_pos(A_new, C_new, S, k1, k2);
      }
      else if (parallel_segments(A_3d, C_3d, S_3d, T_3d)) {
        A_new = xyz2uv(p, A_3d, wall1);
        C_new = xyz2uv(p, C_3d, wall1);
        return get_product_shared_segment_pos(A_new, C_new, R, k1, k2);
      }

    }
    else if (B_shared && C_shared) {
      if (parallel_segments(B_3d, C_3d, R_3d, S_3d)) {
        B_new = xyz2uv(p, B_3d, wall1);
        C_new = xyz2uv(p, C_3d, wall1);
        return get_product_shared_segment_pos(B_new, C_new, T, k1, k2);
      } else if (parallel_segments(B_3d, C_3d, R_3d, T_3d)) {
        B_new = xyz2uv(p, B_3d, wall1);
        C_new = xyz2uv(p, C_3d, wall1);
        return get_product_shared_segment_pos(B_new, C_new, S, k1, k2);
      } else if (parallel_segments(B_3d, C_3d, S_3d, T_3d)) {
        B_new = xyz2uv(p, B_3d, wall1);
        C_new = xyz2uv(p, C_3d, wall1);
        return get_product_shared_segment_pos(B_new, C_new, R, k1, k2);
      }
    }

    /* one vertex shared from each tile */
    if (R_shared) {
      if (A_shared) {
        if (parallel_segments(R_3d, A_3d, R_3d, S_3d)) {
          A_new = xyz2uv(p, A_3d, wall1);
          return get_product_shared_segment_pos(A_new, R, T, k1, k2);
        }
        else if (parallel_segments(R_3d, A_3d, R_3d, T_3d)) {
          A_new = xyz2uv(p, A_3d, wall1);
          return get_product_shared_segment_pos(A_new, R, S, k1, k2);
        }
      }
      else if (B_shared) {
        if (parallel_segments(R_3d, B_3d, R_3d, S_3d)) {
          B_new = xyz2uv(p, B_3d, wall1);
          return get_product_shared_segment_pos(B_new, R, T, k1, k2);
        }
        else if (parallel_segments(R_3d, B_3d, R_3d, T_3d)) {
          B_new = xyz2uv(p, B_3d, wall1);
          return get_product_shared_segment_pos(B_new, R, S, k1, k2);
        }
      }
      else if (C_shared) {
        if (parallel_segments(R_3d, C_3d, R_3d, S_3d)) {
          C_new = xyz2uv(p, C_3d, wall1);
          return get_product_shared_segment_pos(C_new, R, T, k1, k2);
        }
        else if (parallel_segments(R_3d, C_3d, R_3d, T_3d)) {
          C_new = xyz2uv(p, C_3d, wall1);
          return get_product_shared_segment_pos(C_new, R, S, k1, k2);
        }
      }
      else {
        return get_product_shared_vertex_pos(R, S, T, k1, k2);
      }

    }
    else if (S_shared) {
      if (A_shared) {
        if (parallel_segments(S_3d, A_3d, S_3d, T_3d)) {
          A_new = xyz2uv(p, A_3d, wall1);
          return get_product_shared_segment_pos(A_new, S, R, k1, k2);
        }
        else if (parallel_segments(S_3d, A_3d, S_3d, R_3d)) {
          A_new = xyz2uv(p, A_3d, wall1);
          return get_product_shared_segment_pos(A_new, S, T, k1, k2);
        }
      }
      else if (B_shared) {
        if (parallel_segments(S_3d, B_3d, S_3d, T_3d)) {
          B_new = xyz2uv(p, B_3d, wall1);
          return get_product_shared_segment_pos(B_new, S, R, k1, k2);
        }
        else if (parallel_segments(S_3d, B_3d, S_3d, R_3d)) {
          B_new = xyz2uv(p, B_3d, wall1);
          return get_product_shared_segment_pos(B_new, S, T, k1, k2);
        }

      }
      else if (C_shared) {
        if (parallel_segments(S_3d, C_3d, S_3d, T_3d)) {
          C_new = xyz2uv(p, C_3d, wall1);
          return get_product_shared_segment_pos(C_new, S, R, k1, k2);
        }
        else if (parallel_segments(S_3d, C_3d, S_3d, R_3d)) {
          C_new = xyz2uv(p, C_3d, wall1);
          return get_product_shared_segment_pos(C_new, S, T, k1, k2);
        }

      }
      else {
        return get_product_shared_vertex_pos(S, R, T, k1, k2);
      }

    }
    else if (T_shared) {
      if (A_shared) {
        if (parallel_segments(T_3d, A_3d, T_3d, S_3d)) {
          A_new = xyz2uv(p, A_3d, wall1);
          return get_product_shared_segment_pos(A_new, T, R, k1, k2);
        }
        else if (parallel_segments(T_3d, A_3d, T_3d, R_3d)) {
          A_new = xyz2uv(p, A_3d, wall1);
          return get_product_shared_segment_pos(A_new, T, S, k1, k2);
        }
      }
      else if (B_shared) {
        if (parallel_segments(T_3d, B_3d, T_3d, S_3d)) {
          B_new = xyz2uv(p, B_3d, wall1);
          return get_product_shared_segment_pos(B_new, T, R, k1, k2);
        }
        else if (parallel_segments(T_3d, B_3d, T_3d, R_3d)) {
          B_new = xyz2uv(p, B_3d, wall1);
          return get_product_shared_segment_pos(B_new, T, S, k1, k2);
        }

      }
      else if (C_shared) {
        if (parallel_segments(T_3d, C_3d, T_3d, S_3d)) {
          C_new = xyz2uv(p, C_3d, wall1);
          return get_product_shared_segment_pos(C_new, T, R, k1, k2);
        }
        else if (parallel_segments(T_3d, C_3d, T_3d, R_3d)) {
          C_new = xyz2uv(p, C_3d, wall1);
          return get_product_shared_segment_pos(C_new, T, S, k1, k2);
        }

      }
      else {
        return get_product_shared_vertex_pos(T, R, S, k1, k2);
      }
    }

    /* only one vertex is shared */
    if (A_shared) {
      A_new = xyz2uv(p, A_3d, wall1);
      if (intersect_point_segment(A_3d, R_3d, S_3d)) {
        return get_product_close_to_segment_endpoint_pos(T, A_new, k1, k2);
      }
      else if (intersect_point_segment(A_3d, R_3d, T_3d)) {
        return get_product_close_to_segment_endpoint_pos(S, A_new, k1, k2);
      }
      else if (intersect_point_segment(A_3d, S_3d, T_3d)) {
        return get_product_close_to_segment_endpoint_pos(R, A_new, k1, k2);
      }
    }
    else if (B_shared) {
      B_new = xyz2uv(p, B_3d, wall1);
      if (intersect_point_segment(B_3d, R_3d, S_3d)) {
        return get_product_close_to_segment_endpoint_pos(T, B_new, k1, k2);
      }
      else if (intersect_point_segment(B_3d, R_3d, T_3d)) {
        return get_product_close_to_segment_endpoint_pos(S, B_new, k1, k2);
      }
      else if (intersect_point_segment(B_3d, S_3d, T_3d)) {
        return get_product_close_to_segment_endpoint_pos(R, B_new, k1, k2);
      }
    } else if (C_shared) {
      C_new = xyz2uv(p, C_3d, wall1);
      if (intersect_point_segment(C_3d, R_3d, S_3d)) {
        return get_product_close_to_segment_endpoint_pos(T, C_new, k1, k2);
      }
      else if (intersect_point_segment(C_3d, R_3d, T_3d)) {
        return get_product_close_to_segment_endpoint_pos(S, C_new, k1, k2);
      }
      else if (intersect_point_segment(C_3d, S_3d, T_3d)) {
        return get_product_close_to_segment_endpoint_pos(R, C_new, k1, k2);
      }
    }

  } /* end if (num_exact_shared_vertices == 0) */

  /* Apparently there are some round-up errors that force
     the code to come to this place. Below we will try
     again to place the product. */

  /* find points on the triangle RST that are closest to A, B, C */
  Vec3 A_close_3d, B_close_3d, C_close_3d;
  pos_t dist_A_A_close_3d, dist_B_B_close_3d, dist_C_C_close_3d, min_dist;
  Vec3 prod_pos_3d;
  Vec2 prod_pos;

  GeometryUtils::closest_pt_point_triangle(A_3d, R_3d, S_3d, T_3d, A_close_3d);
  GeometryUtils::closest_pt_point_triangle(B_3d, R_3d, S_3d, T_3d, B_close_3d);
  GeometryUtils::closest_pt_point_triangle(C_3d, R_3d, S_3d, T_3d, C_close_3d);

  dist_A_A_close_3d = distance3(A_3d, A_close_3d);
  dist_B_B_close_3d = distance3(B_3d, B_close_3d);
  dist_C_C_close_3d = distance3(C_3d, C_close_3d);

  min_dist = min3_p(dist_A_A_close_3d, dist_B_B_close_3d, dist_C_C_close_3d);

  if (!distinguishable_p(min_dist, dist_A_A_close_3d, POS_EPS)) {
    prod_pos_3d = A_close_3d;
  }
  else if (!distinguishable_p(min_dist, dist_B_B_close_3d, POS_EPS)) {
    prod_pos_3d = B_close_3d;
  }
  else {
    prod_pos_3d = C_close_3d;
  }

  prod_pos = xyz2uv(p, prod_pos_3d, wall1);

  if (intersect_point_segment(prod_pos_3d, R_3d, S_3d)) {
    return get_product_close_to_segment_endpoint_pos(T, prod_pos, k1, k2);
  }
  else if (intersect_point_segment(prod_pos_3d, R_3d, T_3d)) {
    return get_product_close_to_segment_endpoint_pos(S, prod_pos, k1, k2);
  }
  else if (intersect_point_segment(prod_pos_3d, S_3d, T_3d)) {
    return get_product_close_to_segment_endpoint_pos(R, prod_pos, k1, k2);
  }
  else {
    return prod_pos;
  }

  /* I should not come here... */
  mcell_internal_error("Error in the function 'find_closest_position()'.");
}



} // namespace GridPosition
} // namespace MCell

