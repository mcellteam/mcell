#ifndef VECTOR_H
#define VECTOR_H

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
