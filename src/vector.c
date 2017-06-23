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

/* 3D vector routines */

#include "config.h"

#include <math.h>
#include <float.h>

#include "vector.h"
#include "mcell_structs.h"

#define MY_PI 3.14159265358979323846

/**
 * Multiplies two matrices together.
 * Allocates a local 4x4 matrix for computation and then copies result
 * into om.
 * @param m1 a 4 x m input matrix
 * @param m2 a 4 x n input matrix
 * @param om a 4 x n output matrix
 * @param l number of rows in m1
 * @param m number of columns in m1 == number of rows in m2
 * @param n number of columns in m2
 */
void mult_matrix(double (*m1)[4], double (*m2)[4], double (*om)[4],
                 short unsigned int l, short unsigned int m,
                 short unsigned int n) {
  double tm[4][4];
  unsigned short i, j, k;

  for (i = 0; i < l; i++) {
    for (j = 0; j < n; j++) {
      tm[i][j] = 0;
      for (k = 0; k < m; k++) {
        tm[i][j] = tm[i][j] + (m1[i][k]) * (m2[k][j]);
      }
    }
  }
  for (i = 0; i < l; i++) {
    for (j = 0; j < n; j++) {
      om[i][j] = tm[i][j];
    }
  }
}

/**
 * Normalizes a vector3 v.
 */
void normalize(struct vector3 *v) {
  double length;

  length = vect_length(v);
  v->x = v->x / length;
  v->y = v->y / length;
  v->z = v->z / length;
}

/**
 * Initializes a 4x4 Identity matrix.
 */
void init_matrix(double (*im)[4]) {

  im[0][0] = 1;
  im[0][1] = 0;
  im[0][2] = 0;
  im[0][3] = 0;
  im[1][0] = 0;
  im[1][1] = 1;
  im[1][2] = 0;
  im[1][3] = 0;
  im[2][0] = 0;
  im[2][1] = 0;
  im[2][2] = 1;
  im[2][3] = 0;
  im[3][0] = 0;
  im[3][1] = 0;
  im[3][2] = 0;
  im[3][3] = 1;
}

/**
 * Scales the rows of a matrix according to scaling vector3.
 * Scales row0 of im by scale.x
 * Scales row1 of im by scale.y
 * Scales row2 of im by scale.z
 * Scales row3 of im by 1 (no scaling)
 * Result is placed in om.
 */
void scale_matrix(double (*im)[4], double (*om)[4], struct vector3 *scale) {
  double sc[4][4];
  unsigned short l, m, n;

  sc[0][0] = scale->x;
  sc[0][1] = 0;
  sc[0][2] = 0;
  sc[0][3] = 0;
  sc[1][0] = 0;
  sc[1][1] = scale->y;
  sc[1][2] = 0;
  sc[1][3] = 0;
  sc[2][0] = 0;
  sc[2][1] = 0;
  sc[2][2] = scale->z;
  sc[2][3] = 0;
  sc[3][0] = 0;
  sc[3][1] = 0;
  sc[3][2] = 0;
  sc[3][3] = 1;

  l = 4;
  m = 4;
  n = 4;
  mult_matrix(im, sc, om, l, m, n);
}

void translate_matrix(double (*im)[4], double (*om)[4],
                      struct vector3 *translate) {
  double tm[4][4];
  unsigned short l, m, n;

  tm[0][0] = 1;
  tm[0][1] = 0;
  tm[0][2] = 0;
  tm[0][3] = 0;
  tm[1][0] = 0;
  tm[1][1] = 1;
  tm[1][2] = 0;
  tm[1][3] = 0;
  tm[2][0] = 0;
  tm[2][1] = 0;
  tm[2][2] = 1;
  tm[2][3] = 0;
  tm[3][0] = translate->x;
  tm[3][1] = translate->y;
  tm[3][2] = translate->z;
  tm[3][3] = 1;

  l = 4;
  m = 4;
  n = 4;
  mult_matrix(im, tm, om, l, m, n);
}

void rotate_matrix(double (*im)[4], double (*om)[4], struct vector3 *axis,
                   double angle) {
  double r1[4][4], r2[4][4], r3[4][4], rm[4][4];
  double a, b, c, v;
  double rad;
  unsigned short l, m, n;

  normalize(axis);
  a = axis->x;
  b = axis->y;
  c = axis->z;
  v = sqrt(b * b + c * c);

  r1[0][0] = 1;
  r1[0][1] = 0;
  r1[0][2] = 0;
  r1[0][3] = 0;
  r1[1][0] = 0;
  r1[1][1] = 1;
  r1[1][2] = 0;
  r1[1][3] = 0;
  r1[2][0] = 0;
  r1[2][1] = 0;
  r1[2][2] = 1;
  r1[2][3] = 0;
  r1[3][0] = 0;
  r1[3][1] = 0;
  r1[3][2] = 0;
  r1[3][3] = 1;

  if (v != 0.0) {
    r1[1][1] = c / v;
    r1[1][2] = b / v;
    r1[2][1] = -b / v;
    r1[2][2] = c / v;
  }

  r2[0][0] = v;
  r2[0][1] = 0;
  r2[0][2] = a;
  r2[0][3] = 0;
  r2[1][0] = 0;
  r2[1][1] = 1;
  r2[1][2] = 0;
  r2[1][3] = 0;
  r2[2][0] = -a;
  r2[2][1] = 0;
  r2[2][2] = v;
  r2[2][3] = 0;
  r2[3][0] = 0;
  r2[3][1] = 0;
  r2[3][2] = 0;
  r2[3][3] = 1;

  rad = MY_PI / 180.0;
  r3[0][0] = cos(angle * rad);
  r3[0][1] = sin(angle * rad);
  r3[0][2] = 0;
  r3[0][3] = 0;
  r3[1][0] = -sin(angle * rad);
  r3[1][1] = cos(angle * rad);
  r3[1][2] = 0;
  r3[1][3] = 0;
  r3[2][0] = 0;
  r3[2][1] = 0;
  r3[2][2] = 1;
  r3[2][3] = 0;
  r3[3][0] = 0;
  r3[3][1] = 0;
  r3[3][2] = 0;
  r3[3][3] = 1;

  l = 4;
  m = 4;
  n = 4;
  mult_matrix(r1, r2, rm, l, m, n);
  mult_matrix(rm, r3, rm, l, m, n);

  r2[0][2] = -a;
  r2[2][0] = a;

  if (v != 0.0) {
    r1[1][2] = -b / v;
    r1[2][1] = b / v;
  }

  mult_matrix(rm, r2, rm, l, m, n);
  mult_matrix(rm, r1, rm, l, m, n);
  mult_matrix(im, rm, om, l, m, n);
}

void tform_matrix(struct vector3 *scale, struct vector3 *translate,
                  struct vector3 *axis, double angle, double (*om)[4]) {
  double sc[4][4];
  double tm[4][4];
  double r1[4][4], r2[4][4], r3[4][4];
  double a, b, c, v;
  double rad;
  unsigned short l, m, n;

  init_matrix(om);

  sc[0][0] = scale->x;
  sc[0][1] = 0;
  sc[0][2] = 0;
  sc[0][3] = 0;
  sc[1][0] = 0;
  sc[1][1] = scale->y;
  sc[1][2] = 0;
  sc[1][3] = 0;
  sc[2][0] = 0;
  sc[2][1] = 0;
  sc[2][2] = scale->z;
  sc[2][3] = 0;
  sc[3][0] = 0;
  sc[3][1] = 0;
  sc[3][2] = 0;
  sc[3][3] = 1;

  tm[0][0] = 1;
  tm[0][1] = 0;
  tm[0][2] = 0;
  tm[0][3] = 0;
  tm[1][0] = 0;
  tm[1][1] = 1;
  tm[1][2] = 0;
  tm[1][3] = 0;
  tm[2][0] = 0;
  tm[2][1] = 0;
  tm[2][2] = 1;
  tm[2][3] = 0;
  tm[3][0] = translate->x;
  tm[3][1] = translate->y;
  tm[3][2] = translate->z;
  tm[3][3] = 1;

  normalize(axis);
  a = axis->x;
  b = axis->y;
  c = axis->z;
  v = sqrt(b * b + c * c);

  r1[0][0] = 1;
  r1[0][1] = 0;
  r1[0][2] = 0;
  r1[0][3] = 0;
  r1[1][0] = 0;
  r1[1][1] = 1;
  r1[1][2] = 0;
  r1[1][3] = 0;
  r1[2][0] = 0;
  r1[2][1] = 0;
  r1[2][2] = 1;
  r1[2][3] = 0;
  r1[3][0] = 0;
  r1[3][1] = 0;
  r1[3][2] = 0;
  r1[3][3] = 1;

  if (v != 0.0) {
    r1[1][1] = c / v;
    r1[1][2] = b / v;
    r1[2][1] = -b / v;
    r1[2][2] = c / v;
  }

  r2[0][0] = v;
  r2[0][1] = 0;
  r2[0][2] = a;
  r2[0][3] = 0;
  r2[1][0] = 0;
  r2[1][1] = 1;
  r2[1][2] = 0;
  r2[1][3] = 0;
  r2[2][0] = -a;
  r2[2][1] = 0;
  r2[2][2] = v;
  r2[2][3] = 0;
  r2[3][0] = 0;
  r2[3][1] = 0;
  r2[3][2] = 0;
  r2[3][3] = 1;

  rad = MY_PI / 180.0;
  r3[0][0] = cos(angle * rad);
  r3[0][1] = sin(angle * rad);
  r3[0][2] = 0;
  r3[0][3] = 0;
  r3[1][0] = -sin(angle * rad);
  r3[1][1] = cos(angle * rad);
  r3[1][2] = 0;
  r3[1][3] = 0;
  r3[2][0] = 0;
  r3[2][1] = 0;
  r3[2][2] = 1;
  r3[2][3] = 0;
  r3[3][0] = 0;
  r3[3][1] = 0;
  r3[3][2] = 0;
  r3[3][3] = 1;

  l = 4;
  m = 4;
  n = 4;
  mult_matrix(r1, r2, om, l, m, n);
  mult_matrix(om, r3, om, l, m, n);

  r2[0][2] = -a;
  r2[2][0] = a;

  if (v != 0.0) {
    r1[1][2] = -b / v;
    r1[2][1] = b / v;
  }

  mult_matrix(om, r2, om, l, m, n);
  mult_matrix(om, r1, om, l, m, n);
  mult_matrix(om, sc, om, l, m, n);
  mult_matrix(om, tm, om, l, m, n);
}

/**
 * Performs vector subtraction.
 * Subtracts vector3 p1 from vector3 p2 placing the result in vector3 v.
 */
void vectorize(struct vector3 *p1, struct vector3 *p2, struct vector3 *v) {

  v->x = p2->x - p1->x;
  v->y = p2->y - p1->y;
  v->z = p2->z - p1->z;
}

/**
 * Computes the magnitude of a vector.
 */
double vect_length(struct vector3 *v) {
  double length;

  length = sqrt((v->x) * (v->x) + (v->y) * (v->y) + (v->z) * (v->z));
  return (length);
}

/**
 * Computes the dot product of two vector3's v1 and v2.
 */
double dot_prod(struct vector3 *v1, struct vector3 *v2) {
  double dot;

  dot = (v1->x) * (v2->x) + (v1->y) * (v2->y) + (v1->z) * (v2->z);
  return (dot);
}

/**
 * Performs vector cross product.
 * Computes the cross product of two vector3's v1 and v2 storing the result
 * in vector3 v3.
 */
void cross_prod(struct vector3 *v1, struct vector3 *v2, struct vector3 *v3) {

  v3->x = (v1->y) * (v2->z) - (v1->z) * (v2->y);
  v3->y = (v1->z) * (v2->x) - (v1->x) * (v2->z);
  v3->z = (v1->x) * (v2->y) - (v1->y) * (v2->x);
}

/************************************************************************
vect_sum:
 In: v1 - first vector3
     v2 - second vector3
 Out: v3 - the sum of the v1 and v2
************************************************************************/
void vect_sum(struct vector3 *v1, struct vector3 *v2, struct vector3 *v3) {
  v3->x = v1->x + v2->x;
  v3->y = v1->y + v2->y;
  v3->z = v1->z + v2->z;
}
/***********************************************************************
 scalar_prod:
 In: v - vector3
     a - scalar
 Out: result - the product of the vector3 v by scalar a
***********************************************************************/
void scalar_prod(struct vector3 *v1, double a, struct vector3 *result) {
  result->x = a * v1->x;
  result->y = a * v1->y;
  result->z = a * v1->z;
}

/***************************************************************************
distinguishable_vec3 -- reports whether two vectors are measurably different
  (vector analog of distinguishable() in util.c)

Parameters
        a -- first vector
        b -- second vector
        eps -- fractional difference that we think is different

Returns
        1 if the vectors are different, 0 otherwise
***************************************************************************/

int distinguishable_vec3(struct vector3 *a, struct vector3 *b, double eps) {
  double c, cc, d;

  /* Find largest coordinate */
  c = fabs(a->x);

  d = fabs(a->y);
  if (d > c)
    c = d;

  d = fabs(a->z);
  if (d > c)
    c = d;

  d = fabs(b->x);
  if (d > c)
    c = d;

  d = fabs(b->y);
  if (d > c)
    c = d;

  d = fabs(b->z);
  if (d > c)
    c = d;

  /* Find largest difference */
  cc = fabs(a->x - b->x);

  d = fabs(a->y - b->y);
  if (d > cc)
    cc = d;

  d = fabs(a->z - b->z);
  if (d > cc)
    cc = d;

  /* Make sure fractional difference is at least eps and absolute difference is
   * at least (eps*eps) */
  if (c < eps)
    c = eps;
  return (c * eps < cc);
}

/***************************************************************************
distinguishable_vec2 -- reports whether two vectors are measurably different
  (vector analog of distinguishable() in util.c)

Parameters
        a -- first vector2
        b -- second vector2
        eps -- fractional difference that we think is different

Returns
        1 if the vectors are different, 0 otherwise
Note: similar to the function "distinguishable_vec3" but for the surface vectors
***************************************************************************/

int distinguishable_vec2(struct vector2 *a, struct vector2 *b, double eps) {
  double c, cc, d;

  /* Find largest coordinate */
  c = fabs(a->u);

  d = fabs(a->v);
  if (d > c)
    c = d;

  d = fabs(b->u);
  if (d > c)
    c = d;

  d = fabs(b->v);
  if (d > c)
    c = d;

  /* Find largest difference */
  cc = fabs(a->u - b->u);

  d = fabs(a->v - b->v);
  if (d > cc)
    cc = d;

  /* Make sure fractional difference is at least eps and absolute difference is
   * at least (eps*eps) */
  if (c < eps)
    c = eps;
  return (c * eps < cc);
}

/***************************************************************************
distance_vec3 -- calculates distance between two points in 3D

Parameters
        a -- first point
        b -- second point

Returns
        distance between two points in 3D
***************************************************************************/
double distance_vec3(struct vector3 *a, struct vector3 *b) {
  double dist;
  dist = sqrt((a->x - b->x) * (a->x - b->x) + (a->y - b->y) * (a->y - b->y) +
              (a->z - b->z) * (a->z - b->z));

  return dist;
}

/****************************************************************************
parallel_segments:
   In: segment defined by endpoints A, B
       segment defined by endpoints R, S
   Out: 1, if the segments are parallel.
        0, otherwise
*****************************************************************************/
int parallel_segments(struct vector3 *A, struct vector3 *B, struct vector3 *R,
                      struct vector3 *S) {

  double length;
  struct vector3 prod; /* cross product */
  struct vector3 ba, sr;

  vectorize(A, B, &ba);
  vectorize(S, R, &sr);
  cross_prod(&ba, &sr, &prod);

  length = vect_length(&prod);

  if (!distinguishable(length, 0, EPS_C))
    return 1;

  return 0;
}
/**************************************************************************
same_side:
        In: two points p1 and p2
            line defined by the points a and b
        Out: returns 1 if points p1 and p2 are on the same side of the line
             defined by the points a and b
**************************************************************************/
int same_side(struct vector3 *p1, struct vector3 *p2, struct vector3 *a,
              struct vector3 *b) {
  struct vector3 cp1, cp2, b_a, p1_a, p2_a;
  vectorize(a, b, &b_a);
  vectorize(a, p1, &p1_a);
  vectorize(a, p2, &p2_a);
  cross_prod(&b_a, &p1_a, &cp1);
  cross_prod(&b_a, &p2_a, &cp2);

  if (dot_prod(&cp1, &cp2) >= 0) {
    return 1;
  } else
    return 0;
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
int point_in_triangle(struct vector3 *p, struct vector3 *a, struct vector3 *b,
                      struct vector3 *c) {
  if (same_side(p, a, b, c) && same_side(p, b, a, c) && same_side(p, c, a, b)) {
    return 1;
  }

  if (((!distinguishable(p->x, a->x, EPS_C)) &&
       (!distinguishable(p->y, a->y, EPS_C)) &&
       (!distinguishable(p->z, a->z, EPS_C))) ||
      ((!distinguishable(p->x, b->x, EPS_C)) &&
       (!distinguishable(p->y, b->y, EPS_C)) &&
       (!distinguishable(p->z, b->z, EPS_C))) ||
      ((!distinguishable(p->x, c->x, EPS_C)) &&
       (!distinguishable(p->y, c->y, EPS_C)) &&
       (!distinguishable(p->z, c->z, EPS_C)))) {
    return 1;
  }

  return 0;
}

#undef MY_PI

/*******************************************************************
cross2D:
   In: 2D vectors a and b
   Out: 2D pseudo cross product Dot(Perp(a0,b)
   Note: The code adapted from "Real-Time Collision Detection" by
              Christer Ericson, p.205

*******************************************************************/
double cross2D(struct vector2 *a, struct vector2 *b) {
  return ((a->v) * (b->u) - (a->u) * (b->v));
}

/*************************************************************************
vectorize2D:
   In: 2D vectors p1 and p2
   Out: Subtracts vector p1 from p2 and places result into p3
*************************************************************************/
void vectorize2D(struct vector2 *p1, struct vector2 *p2, struct vector2 *p3) {
  p3->u = p2->u - p1->u;
  p3->v = p2->v - p1->v;
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
int point_in_triangle_2D(struct vector2 *p, struct vector2 *a,
                         struct vector2 *b, struct vector2 *c) {
  struct vector2 p_minus_a, b_minus_a, p_minus_b, c_minus_b, p_minus_c,
      a_minus_c;
  double pab, pbc, pca;

  vectorize2D(a, p, &p_minus_a);
  vectorize2D(a, b, &b_minus_a);
  vectorize2D(b, p, &p_minus_b);
  vectorize2D(b, c, &c_minus_b);
  vectorize2D(c, p, &p_minus_c);
  vectorize2D(c, a, &a_minus_c);

  pab = cross2D(&p_minus_a, &b_minus_a);
  pbc = cross2D(&p_minus_b, &c_minus_b);
  /* if P left of one of AB and BC and right of the other, not inside triangle
     - (pab and pbc have different signs */
  if (((pab > 0) && (pbc < 0)) || ((pab < 0) && (pbc > 0)))
    return 0;

  pca = cross2D(&p_minus_c, &a_minus_c);
  /* if P left of one of AB and CA and right of the other, not inside triangle
   * - pab and pca have different signs */
  if (((pab > 0) && (pca < 0)) || ((pab < 0) && (pca > 0)))
    return 0;

  /* if P left or right of all edges, so must be in (or on) the triangle */
  return 1;
}

/*****************************************************************************
intersect_point_segment:
   In: point P and segment AB
   Out: 1 if the point P lies on the segment AB, and 0 - otherwise
******************************************************************************/
int intersect_point_segment(struct vector3 *P, struct vector3 *A,
                            struct vector3 *B) {
  struct vector3 ba, pa;
  double t;                    /* parameter in the line parametrical equation */
  double ba_length, pa_length; /* length of the vectors */
  double cosine_angle;         /* cosine of the angle between ba and pa */

  /* check for the end points */
  if (!distinguishable_vec3(P, A, EPS_C))
    return 1;
  if (!distinguishable_vec3(P, B, EPS_C))
    return 1;

  vectorize(A, B, &ba);
  vectorize(A, P, &pa);

  ba_length = vect_length(&ba);
  pa_length = vect_length(&pa);

  /* if point intersects segment, vectors pa and ba should be collinear */
  cosine_angle = dot_prod(&ba, &pa) / (ba_length * pa_length);
  if (distinguishable(cosine_angle, 1.0, EPS_C)) {
    return 0;
  }

  /* Project P on AB, computing parameterized position d(t) = A + t(B - A) */
  t = dot_prod(&pa, &ba) / dot_prod(&ba, &ba);

  if ((t > 0) && (t < 1))
    return 1;

  return 0;
}

/***************************************************************************
point_in_box:
  In:  low_left - lower left front corner of the box
       up_right - upper right back corner of the box
       point - we want to find out  if this point is in the box
  Out: Returns 1 if point is in box, 0 otherwise
       This is very similar to test_bounding_boxes.
***************************************************************************/
int point_in_box(struct vector3 *low_left, struct vector3 *up_right,
                 struct vector3 *point) {
  if ((up_right->x < point->x) || (low_left->x > point->x))
    return 0;
  if ((up_right->y < point->y) || (low_left->y > point->y))
    return 0;
  if ((up_right->z < point->z) || (low_left->z > point->z))
    return 0;
  return 1;
}
