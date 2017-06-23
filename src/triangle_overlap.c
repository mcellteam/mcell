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

/* The triangle overlap code in this file is based on the triangle/triangle
 * intersection test routine by by Tomas Moller, 1997.
 *
 * See article "A Fast Triangle-Triangle Intersection Test",
 * Journal of Graphics Tools, 2(2), 1997
 *
 * In contrast to Moller's general triangle-triangle intersection routine
 * our code only tests for area overlap of coplanar triangles.
 * In particular it will treat triangles which share an edge or vertex as
 * non-overlapping. Only triangles with a bona-fide area overlap will be
 * flagged. 
 *
 * The code returns 1 if the provided triangles have an area overlap and 0
 * otherwise.
 */

#include <stdbool.h>
#include <math.h>

#include "triangle_overlap.h"

#define EPSILON EPS_C

/* some macros */
#define CROSS(dest, v1, v2)                  \
  {                                          \
    dest[0] = v1[1] * v2[2] - v1[2] * v2[1]; \
    dest[1] = v1[2] * v2[0] - v1[0] * v2[2]; \
    dest[2] = v1[0] * v2[1] - v1[1] * v2[0]; \
  }

#define SUB(dest, v1, v2)    \
  {                          \
    dest[0] = v1[0] - v2[0]; \
    dest[1] = v1[1] - v2[1]; \
    dest[2] = v1[2] - v2[2]; \
  }

/* this edge to edge test is based on Franlin Antonio's gem:
 * "Faster Line Segment Intersection", in Graphics Gems III,
 * pp. 199-202 */
static inline int edge_edge_test(double* V0, double* U0, double* U1, short i0,
                                 short i1, double Ax, double Ay) {
  double Bx = U0[i0] - U1[i0];
  double By = U0[i1] - U1[i1];
  double Cx = V0[i0] - U0[i0];
  double Cy = V0[i1] - U0[i1];
  double f = Ay * Bx - Ax * By;
  double d = By * Cx - Bx * Cy;
  if ((f > 0 && d >= 0 && d <= f) || (f < 0 && d <= 0 && d >= f)) {
    double e = Ax * Cy - Ay * Cx;

    // ignore edge or vertex overlaps
    bool dz = !distinguishable(d, 0.0, EPSILON);
    bool ez = !distinguishable(e, 0.0, EPSILON);
    bool df = !distinguishable(d, f, EPSILON);
    bool ef = !distinguishable(e, f, EPSILON);
    if ((dz && ez) || (dz && ef) || (df && ez) || (df && ef)) {
      return 0;
    }

    if (f > 0) {
      if (e >= 0 && e <= f) return 1;
    } else {
      if (e <= 0 && e >= f) return 1;
    }
  }
  return 0;
}

/* edge_against_tri_edges tests if edge V0V1 intersects with any edge of 
 * triangle U0U1U2 */
static inline int edge_against_tri_edges(double* V0, double* V1, double* U0,
                                         double* U1, double* U2, short i0,
                                         short i1) {
  double Ax = V1[i0] - V0[i0];
  double Ay = V1[i1] - V0[i1];
  /* test edge U0,U1 against V0,V1 */
  if (edge_edge_test(V0, U0, U1, i0, i1, Ax, Ay)) {
    return 1;
  }
  /* test edge U1,U2 against V0,V1 */
  if (edge_edge_test(V0, U1, U2, i0, i1, Ax, Ay)) {
    return 1;
  };
  /* test edge U2,U1 against V0,V1 */
  if (edge_edge_test(V0, U2, U0, i0, i1, Ax, Ay)) {
    return 1;
  }
  return 0;
}

/* point_in_tri tests if point V0 is contained in triangle U0U1U2 */
static inline int point_in_tri(double* V0, double* U0, double* U1, double* U2,
                               short i0, short i1) {
  /* is T1 completly inside T2? */
  /* check if V0 is inside tri(U0,U1,U2) */
  double a = U1[i1] - U0[i1];
  double b = -(U1[i0] - U0[i0]);
  double c = -a * U0[i0] - b * U0[i1];
  double d0 = a * V0[i0] + b * V0[i1] + c;

  a = U2[i1] - U1[i1];
  b = -(U2[i0] - U1[i0]);
  c = -a * U1[i0] - b * U1[i1];
  double d1 = a * V0[i0] + b * V0[i1] + c;

  a = U0[i1] - U2[i1];
  b = -(U0[i0] - U2[i0]);
  c = -a * U2[i0] - b * U2[i1];
  double d2 = a * V0[i0] + b * V0[i1] + c;
  if (d0 * d1 > 0.0) {
    if (d0 * d2 > 0.0) return 1;
  }
  return 0;
}

/* coplanar_tri_tri checks if triangles V0V1V2 and U0U1U2 have area overlap. */
static inline int coplanar_tri_tri(double* N, double* V0, double* V1,
                                   double* V2, double* U0, double* U1,
                                   double* U2) {
  double A[3];
  short i0, i1;
  /* first project onto an axis-aligned plane, that maximizes the area */
  /* of the triangles, compute indices: i0,i1. */
  A[0] = fabs(N[0]);
  A[1] = fabs(N[1]);
  A[2] = fabs(N[2]);
  if (A[0] > A[1]) {
    if (A[0] > A[2]) {
      i0 = 1; /* A[0] is greatest */
      i1 = 2;
    } else {
      i0 = 0; /* A[2] is greatest */
      i1 = 1;
    }
  } else /* A[0]<=A[1] */
  {
    if (A[2] > A[1]) {
      i0 = 0; /* A[2] is greatest */
      i1 = 1;
    } else {
      i0 = 0; /* A[1] is greatest */
      i1 = 2;
    }
  }

  /* test all edges of triangle 1 against the edges of triangle 2 */
  if (edge_against_tri_edges(V0, V1, U0, U1, U2, i0, i1) ||
      edge_against_tri_edges(V1, V2, U0, U1, U2, i0, i1) ||
      edge_against_tri_edges(V2, V0, U0, U1, U2, i0, i1)) {
    return 1;
  }

  /* finally, test if tri1 is totally contained in tri2 or vice versa */
  if (point_in_tri(V0, U0, U1, U2, i0, i1) ||
      point_in_tri(U0, V0, V1, V2, i0, i1)) {
    return 1;
  }
  return 0;
}

/* coplanar_tri_operlap tests if w1 and w2 overlap and returns 1 if yes
 * and 0 otherwise.
 *
 * NOTE 1: w1 and w2 are assumed to be coplanar.
 * NOTE 2: this function will only return 1 if w1 and w2 have a bona
 *         fide area overlap. If w1 and w2 share an edge or vertex
 *         this will *not* be treated as an overlap.
 */
int coplanar_tri_overlap(struct wall* w1, struct wall* w2) {
  double V0[3] = {w1->vert[0]->x, w1->vert[0]->y, w1->vert[0]->z};
  double V1[3] = {w1->vert[1]->x, w1->vert[1]->y, w1->vert[1]->z};
  double V2[3] = {w1->vert[2]->x, w1->vert[2]->y, w1->vert[2]->z};
  double U0[3] = {w2->vert[0]->x, w2->vert[0]->y, w2->vert[0]->z};
  double U1[3] = {w2->vert[1]->x, w2->vert[1]->y, w2->vert[1]->z};
  double U2[3] = {w2->vert[2]->x, w2->vert[2]->y, w2->vert[2]->z};

  /* compute plane equation of triangle(V0,V1,V2) */
  double E1[3], E2[3];
  SUB(E1, V1, V0);
  SUB(E2, V2, V0);

  double N1[3];
  CROSS(N1, E1, E2);
  /* plane equation 1: N1.X+d1=0 */
  return coplanar_tri_tri(N1, V0, V1, V2, U0, U1, U2);
}
