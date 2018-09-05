/******************************************************************************
 *
 * Copyright (C) 2006-2015 by
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


/* Header file for 3D vector routines */

struct vector2 {
  double u;
  double v;
};

struct vector3 {
  double x;
  double y;
  double z;
};

void mult_matrix(double (*m1)[4], double (*m2)[4], double (*om)[4],
                 short unsigned int l, short unsigned int m,
                 short unsigned int n);
void normalize(struct vector3 *v);
void init_matrix(double (*im)[4]);
void scale_matrix(double (*im)[4], double (*om)[4], struct vector3 *scale);
void translate_matrix(double (*im)[4], double (*om)[4],
                      struct vector3 *translate);
void rotate_matrix(double (*im)[4], double (*om)[4], struct vector3 *axis,
                   double angle);
void tform_matrix(struct vector3 *scale, struct vector3 *translate,
                  struct vector3 *axis, double angle, double (*om)[4]);
void vectorize(struct vector3 *p1, struct vector3 *p2, struct vector3 *v);
double vect_length(struct vector3 *v);
double dot_prod(struct vector3 *v1, struct vector3 *v2);
void cross_prod(struct vector3 *v1, struct vector3 *v2, struct vector3 *v3);
void vect_sum(struct vector3 *v1, struct vector3 *v2, struct vector3 *v3);
void scalar_prod(struct vector3 *v1, double a, struct vector3 *result);

int distinguishable_vec3(struct vector3 *a, struct vector3 *b, double eps);
int distinguishable_vec2(struct vector2 *a, struct vector2 *b, double eps);

double distance_vec3(struct vector3 *a, struct vector3 *b);

int parallel_segments(struct vector3 *A, struct vector3 *B, struct vector3 *R,
                      struct vector3 *S);

int point_in_triangle(struct vector3 *p, struct vector3 *a, struct vector3 *b,
                      struct vector3 *c);
int same_side(struct vector3 *p1, struct vector3 *p2, struct vector3 *a,
              struct vector3 *b);

int intersect_point_segment(struct vector3 *P, struct vector3 *A,
                            struct vector3 *B);

double cross2D(struct vector2 *a, struct vector2 *b);
void vectorize2D(struct vector2 *p1, struct vector2 *p2, struct vector2 *p3);
int point_in_triangle_2D(struct vector2 *p, struct vector2 *a,
                         struct vector2 *b, struct vector2 *c);

int point_in_box(struct vector3 *low_left, struct vector3 *up_right,
                 struct vector3 *point);
