#ifndef VECTOR_H
#define VECTOR_H

/*
MCell (tm) Version 2.08 5/18/1998

Copyright (C) 1997,1998, The Salk Institute & Cornell University.
MCell was written jointly by T.M. Bartol Jr. & J.R. Stiles,
with design input for Monte Carlo algorithms from E.E. Salpeter.

  Acknowledgements:
    T.J. Sejnowski for development input and support
    (NSF Grant IBN-9603611), and M.M. Salpeter for fostering
    quantitative experimental applications.  Additional support
    from NIH Grant K08NS01776 (J.R. Stiles).

MCell is a scientific simulation software tool distributed freely to
registered beta-test research laboratories, and must be obtained as a
machine-specific executable program from the authors.  Copying and/or
editing of MCell without the authors' consent is expressly forbidden.
MCell is provided as is and the authors assume no responsibility for
simulation results obtained by users.  Any material published from
MCell simulations must acknowledge the authors and granting
institutions.
*/

/* Header file for 3D vector routines */

#ifndef VECTOR_MATH
#define VECTOR_MATH

struct vector2 {
	double u;
	double v;
};

struct vector3 {
	double x;
	double y;
	double z;
};

void mult_matrix(double (*m1)[4], double (*m2)[4], double (*om)[4], short unsigned int l, short unsigned int m, short unsigned int n);
void normalize(struct vector3 *v);
void init_matrix(double (*im)[4]);
void scale_matrix(double (*im)[4], double (*om)[4], struct vector3 *scale);
void translate_matrix(double (*im)[4], double (*om)[4], struct vector3 *translate);
void rotate_matrix(double (*im)[4], double (*om)[4], struct vector3 *axis, double angle);
void tform_matrix(struct vector3 *scale, struct vector3 *translate, struct vector3 *axis, double angle, double (*om)[4]);
void vectorize(struct vector3 *p1, struct vector3 *p2, struct vector3 *v);
double vect_length(struct vector3 *v);
double dot_prod(struct vector3 *v1, struct vector3 *v2);
void cross_prod(struct vector3 *v1, struct vector3 *v2, struct vector3 *v3);

#endif

#endif
