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

/* 3D vector routines */

#include <math.h>
#include "vector.h"

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
void mult_matrix(double (*m1)[4], double (*m2)[4], double (*om)[4], short unsigned int l, short unsigned int m, short unsigned int n)
{
  double tm[4][4];
  unsigned short i,j,k;

  for (i=0;i<l;i++) {
    for (j=0;j<n;j++) {
      tm[i][j]=0;
      for (k=0;k<m;k++) {
	tm[i][j]=tm[i][j]+(m1[i][k])*(m2[k][j]);
      }
    }
  }
  for (i=0;i<l;i++) {
    for (j=0;j<n;j++) {
      om[i][j]=tm[i][j];
    }
  }
}


/**
 * Normalizes a vector3 v.
 */
void normalize(struct vector3 *v)
{
  double length;

  length=vect_length(v);
  v->x=v->x/length;
  v->y=v->y/length;
  v->z=v->z/length;
}


/**
 * Initializes a 4x4 Identity matrix.
 */
void init_matrix(double (*im)[4])
{

  im[0][0]=1;
  im[0][1]=0;
  im[0][2]=0;
  im[0][3]=0;
  im[1][0]=0;
  im[1][1]=1;
  im[1][2]=0;
  im[1][3]=0;
  im[2][0]=0;
  im[2][1]=0;
  im[2][2]=1;
  im[2][3]=0;
  im[3][0]=0;
  im[3][1]=0;
  im[3][2]=0;
  im[3][3]=1;
}

/**
 * Scales the rows of a matrix according to scaling vector3.
 * Scales row0 of im by scale.x
 * Scales row1 of im by scale.y
 * Scales row2 of im by scale.z
 * Scales row3 of im by 1 (no scaling)
 * Result is placed in om.
 */
void scale_matrix(double (*im)[4], double (*om)[4], struct vector3 *scale)
{
  double sc[4][4];
  unsigned short l,m,n;

  sc[0][0]=scale->x;
  sc[0][1]=0;
  sc[0][2]=0;
  sc[0][3]=0;
  sc[1][0]=0;
  sc[1][1]=scale->y;
  sc[1][2]=0;
  sc[1][3]=0;
  sc[2][0]=0;
  sc[2][1]=0;
  sc[2][2]=scale->z;
  sc[2][3]=0;
  sc[3][0]=0;
  sc[3][1]=0;
  sc[3][2]=0;
  sc[3][3]=1;

  l=4;
  m=4;
  n=4;
  mult_matrix(im,sc,om,l,m,n);
}

void translate_matrix(double (*im)[4], double (*om)[4], struct vector3 *translate)
{
  double tm[4][4];
  unsigned short l,m,n;

  tm[0][0]=1;
  tm[0][1]=0;
  tm[0][2]=0;
  tm[0][3]=0;
  tm[1][0]=0;
  tm[1][1]=1;
  tm[1][2]=0;
  tm[1][3]=0;
  tm[2][0]=0;
  tm[2][1]=0;
  tm[2][2]=1;
  tm[2][3]=0;
  tm[3][0]=translate->x;
  tm[3][1]=translate->y;
  tm[3][2]=translate->z;
  tm[3][3]=1;

  l=4;
  m=4;
  n=4;
  mult_matrix(im,tm,om,l,m,n);
}

void rotate_matrix(double (*im)[4], double (*om)[4], struct vector3 *axis, double angle)
{
  double r1[4][4],r2[4][4],r3[4][4],rm[4][4];
  double a,b,c,v;
  double rad;
  unsigned short l,m,n;

  normalize(axis);
  a=axis->x;
  b=axis->y;
  c=axis->z;
  v=sqrt(b*b+c*c);

  r1[0][0]=1;
  r1[0][1]=0;
  r1[0][2]=0;
  r1[0][3]=0;
  r1[1][0]=0;
  r1[1][1]=1;
  r1[1][2]=0;
  r1[1][3]=0;
  r1[2][0]=0;
  r1[2][1]=0;
  r1[2][2]=1;
  r1[2][3]=0;
  r1[3][0]=0;
  r1[3][1]=0;
  r1[3][2]=0;
  r1[3][3]=1;

  if (v!=0.0) {
    r1[1][1]=c/v;
    r1[1][2]=b/v;
    r1[2][1]=-b/v;
    r1[2][2]=c/v;
  }

  r2[0][0]=v;
  r2[0][1]=0;
  r2[0][2]=a;
  r2[0][3]=0;
  r2[1][0]=0;
  r2[1][1]=1;
  r2[1][2]=0;
  r2[1][3]=0;
  r2[2][0]=-a;
  r2[2][1]=0;
  r2[2][2]=v;
  r2[2][3]=0;
  r2[3][0]=0;
  r2[3][1]=0;
  r2[3][2]=0;
  r2[3][3]=1;

  rad=MY_PI/180.0;
  r3[0][0]=cos(angle*rad);
  r3[0][1]=sin(angle*rad);
  r3[0][2]=0;
  r3[0][3]=0;
  r3[1][0]=-sin(angle*rad);
  r3[1][1]=cos(angle*rad);
  r3[1][2]=0;
  r3[1][3]=0;
  r3[2][0]=0;
  r3[2][1]=0;
  r3[2][2]=1;
  r3[2][3]=0;
  r3[3][0]=0;
  r3[3][1]=0;
  r3[3][2]=0;
  r3[3][3]=1;

  l=4;
  m=4;
  n=4;
  mult_matrix(r1,r2,rm,l,m,n);
  mult_matrix(rm,r3,rm,l,m,n);

  r2[0][2]=-a;
  r2[2][0]=a;

  if (v!=0.0) {
    r1[1][2]=-b/v;
    r1[2][1]=b/v;
  }

  mult_matrix(rm,r2,rm,l,m,n);
  mult_matrix(rm,r1,rm,l,m,n);
  mult_matrix(im,rm,om,l,m,n);
}

void tform_matrix(struct vector3 *scale, struct vector3 *translate, struct vector3 *axis, double angle, double (*om)[4])
{
  double sc[4][4];
  double tm[4][4];
  double r1[4][4],r2[4][4],r3[4][4];
  double a,b,c,v;
  double rad;
  unsigned short l,m,n;

  init_matrix(om);

  sc[0][0]=scale->x;
  sc[0][1]=0;
  sc[0][2]=0;
  sc[0][3]=0;
  sc[1][0]=0;
  sc[1][1]=scale->y;
  sc[1][2]=0;
  sc[1][3]=0;
  sc[2][0]=0;
  sc[2][1]=0;
  sc[2][2]=scale->z;
  sc[2][3]=0;
  sc[3][0]=0;
  sc[3][1]=0;
  sc[3][2]=0;
  sc[3][3]=1;

  tm[0][0]=1;
  tm[0][1]=0;
  tm[0][2]=0;
  tm[0][3]=0;
  tm[1][0]=0;
  tm[1][1]=1;
  tm[1][2]=0;
  tm[1][3]=0;
  tm[2][0]=0;
  tm[2][1]=0;
  tm[2][2]=1;
  tm[2][3]=0;
  tm[3][0]=translate->x;
  tm[3][1]=translate->y;
  tm[3][2]=translate->z;
  tm[3][3]=1;

  normalize(axis);
  a=axis->x;
  b=axis->y;
  c=axis->z;
  v=sqrt(b*b+c*c);

  r1[0][0]=1;
  r1[0][1]=0;
  r1[0][2]=0;
  r1[0][3]=0;
  r1[1][0]=0;
  r1[1][1]=1;
  r1[1][2]=0;
  r1[1][3]=0;
  r1[2][0]=0;
  r1[2][1]=0;
  r1[2][2]=1;
  r1[2][3]=0;
  r1[3][0]=0;
  r1[3][1]=0;
  r1[3][2]=0;
  r1[3][3]=1;

  if (v!=0.0) {
    r1[1][1]=c/v;
    r1[1][2]=b/v;
    r1[2][1]=-b/v;
    r1[2][2]=c/v;
  }

  r2[0][0]=v;
  r2[0][1]=0;
  r2[0][2]=a;
  r2[0][3]=0;
  r2[1][0]=0;
  r2[1][1]=1;
  r2[1][2]=0;
  r2[1][3]=0;
  r2[2][0]=-a;
  r2[2][1]=0;
  r2[2][2]=v;
  r2[2][3]=0;
  r2[3][0]=0;
  r2[3][1]=0;
  r2[3][2]=0;
  r2[3][3]=1;

  rad=MY_PI/180.0;
  r3[0][0]=cos(angle*rad);
  r3[0][1]=sin(angle*rad);
  r3[0][2]=0;
  r3[0][3]=0;
  r3[1][0]=-sin(angle*rad);
  r3[1][1]=cos(angle*rad);
  r3[1][2]=0;
  r3[1][3]=0;
  r3[2][0]=0;
  r3[2][1]=0;
  r3[2][2]=1;
  r3[2][3]=0;
  r3[3][0]=0;
  r3[3][1]=0;
  r3[3][2]=0;
  r3[3][3]=1;

  l=4;
  m=4;
  n=4;
  mult_matrix(r1,r2,om,l,m,n);
  mult_matrix(om,r3,om,l,m,n);

  r2[0][2]=-a;
  r2[2][0]=a;

  if (v!=0.0) {
    r1[1][2]=-b/v;
    r1[2][1]=b/v;
  }

  mult_matrix(om,r2,om,l,m,n);
  mult_matrix(om,r1,om,l,m,n);
  mult_matrix(om,sc,om,l,m,n);
  mult_matrix(om,tm,om,l,m,n);
}


/**
 * Performs vector subtraction.
 * Subtracts vector3 p1 from vector3 p2 placing the result in vector3 v.
 */
void vectorize(struct vector3 *p1, struct vector3 *p2, struct vector3 *v)
{

  v->x=p2->x-p1->x;
  v->y=p2->y-p1->y;
  v->z=p2->z-p1->z;
}


/**
 * Computes the magnitude of a vector.
 */
double vect_length(struct vector3 *v)
{
  double length;

  length=sqrt((v->x)*(v->x)+(v->y)*(v->y)+(v->z)*(v->z));
  return(length);
}

/**
 * Computes the dot product of two vector3's v1 and v2.
 */
double dot_prod(struct vector3 *v1, struct vector3 *v2)
{
  double dot;

  dot=(v1->x)*(v2->x)+(v1->y)*(v2->y)+(v1->z)*(v2->z);
  return(dot);
}

/**
 * Performs vector cross product.
 * Computes the cross product of two vector3's v1 and v2 storing the result
 * in vector3 v3.
 */
void cross_prod(struct vector3 *v1, struct vector3 *v2, struct vector3 *v3)
{

  v3->x=(v1->y)*(v2->z)-(v1->z)*(v2->y);
  v3->y=(v1->z)*(v2->x)-(v1->x)*(v2->z);
  v3->z=(v1->x)*(v2->y)-(v1->y)*(v2->x);
}

/************************************************************************
vect_sum:
 In: v1 - first vector3
     v2 - second vector3
 Out: v3 - the sum of the v1 and v2 
************************************************************************/
void vect_sum(struct vector3 *v1, struct vector3 *v2, struct vector3 *v3)
{
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
void scalar_prod(struct vector3 *v1, double a, struct vector3 *result)
{
  result->x = a*v1->x;
  result->y = a*v1->y;
  result->z = a*v1->z;
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

int distinguishable_vec3(struct vector3 *a,struct vector3 *b,double eps)
{
  double c,cc,d;
  
  /* Find largest coordinate */
  c=a->x;
  if (c<0) c=-c;
  
  d=a->y;
  if (d<0) d=-d;
  if (d>c) c=d;
  
  d=a->z;
  if (d<0) d=-d;
  if (d>c) c=d;

  d=b->x;
  if (d<0) d=-d;
  if (d>c) c=d;
  
  d=b->y;
  if (d<0) d=-d;
  if (d>c) c=d;
  
  d=b->z;
  if (d<0) d=-d;
  if (d>c) c=d;
  
  /* Find largest difference */
  cc=a->x-b->x;
  if (cc<0) cc=-cc;
  
  d=a->y-b->y;
  if (d<0) d=-d;
  if (d>cc) cc=d;
  
  d=a->z-b->z;
  if (d<0) d=-d;
  if (d>cc) cc=d;
  
  /* Make sure fractional difference is at least eps and absolute difference is at least (eps*eps) */
  if (d<eps) d=eps;
  return (d*eps < c);
}

#undef MY_PI
