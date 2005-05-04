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

/*
 * Adapted from Numerical Recipes in C
 * (C) Copr. 1986-92 Numerical Recipes Software #.,. 
 */
#include "ran4.h"
#include <math.h>
#include <stdlib.h>

#define NITER 2


/**
 * Initialize internal state of random number generator
 * with seed value coming from external state variable.
 * Initialize external state variable to 1.
 */
void ran4_init(struct ran4_state *rng,int seed)
{
  rng->i_state = seed;
  rng->e_state = 1;
}


/**
 * Generate an array of double precision random numbers.
 */
double ran4_dbl32(struct ran4_state *rng)
{
  register unsigned int irword,lword;
  unsigned int j;
  register unsigned int ia,ib,iswap,itmph=0,itmpl=0;
  static unsigned int c1[4]={
	  0xbaa96887L, 0x1e17d32cL, 0x03bcdc3cL, 0x0f33d1b2L};
  static unsigned int c2[4]={
	  0x4b0f3b58L, 0xe874f0c3L, 0x6955c5a6L, 0x55a7ca46L};

  static double range=2.3283064365386963e-10;
  lword=rng->i_state;
  irword=rng->e_state++;

  for (j=0;j<NITER;j++) {
	ia=(iswap=(irword)) ^ c1[j];
	itmpl = ia & 0xffff;
	itmph = ia >> 16;
	ib=itmpl*itmpl+ ~(itmph*itmph);
	irword=(lword) ^ (((ia = (ib >> 16) |
		((ib & 0xffff) << 16)) ^ c2[j])+itmpl*itmph);
	lword=iswap;
  }
  
  return range*irword;
}



/**
 * Generate an array of integer random numbers.
 */
unsigned int ran4_uint32(struct ran4_state *rng)
{
  unsigned int irword,lword;
  unsigned int j,ia,ib,iswap,itmph=0,itmpl=0;
  static unsigned int c1[4]={
	  0xbaa96887L, 0x1e17d32cL, 0x03bcdc3cL, 0x0f33d1b2L};
  static unsigned int c2[4]={
	  0x4b0f3b58L, 0xe874f0c3L, 0x6955c5a6L, 0x55a7ca46L};

  lword=rng->i_state;
  irword=rng->e_state++;

  for (j=0;j<NITER;j++) {
	ia=(iswap=(irword)) ^ c1[j];
	itmpl = ia & 0xffff;
	itmph = ia >> 16;
	ib=itmpl*itmpl+ ~(itmph*itmph);
	irword=(lword) ^ (((ia = (ib >> 16) |
		((ib & 0xffff) << 16)) ^ c2[j])+itmpl*itmph);
	lword=iswap;
  }

  return irword;
}


#if 0

double ran4_dbl32_gauss(struct ran4_state *rng)
{
  double tmp_vec[2];
  double fac,r,v1,v2;
  unsigned int i;

  v1 = ran4_dbl32(rng) - 1.0;
  v2 = ran4_dbl32(rng) - 1.0;
  r = v1*v1 + v2*v2;
  
  return v1*sqrt(-2.0*log(r)/r);
}


/* Generate a single random unsigned integer using ran4 method */

unsigned int rng_uint(unsigned int index)
{
  unsigned int lword;
  unsigned int j,ia,ib,iswap,itmph=0,itmpl=0;
  static unsigned int c1[4]={0xbaa96887L, 0x1e17d32cL, 0x03bcdc3cL, 0x0f33d1b2L};
  static unsigned int c2[4]={0x4b0f3b58L, 0xe874f0c3L, 0x6955c5a6L, 0x55a7ca46L};

  lword=i_state;

  for (j=0;j<NITER;j++)
  {
    ia=(iswap=(index)) ^ c1[j];
    itmpl = ia & 0xffff;
    itmph = ia >> 16;
    ib=itmpl*itmpl+ ~(itmph*itmph);
    index=(lword) ^ (((ia = (ib >> 16) |
                      ((ib & 0xffff) << 16)) ^ c2[j])+itmpl*itmph);
    lword=iswap;
  }
  
  return index;
}


/* Generate a single random double using ran4 method */

double rng_double(unsigned int index)
{
  unsigned int lword;
  unsigned int j,ia,ib,iswap,itmph=0,itmpl=0;
  static unsigned int c1[4]={0xbaa96887L, 0x1e17d32cL, 0x03bcdc3cL, 0x0f33d1b2L};
  static unsigned int c2[4]={0x4b0f3b58L, 0xe874f0c3L, 0x6955c5a6L, 0x55a7ca46L};

  lword=i_state;

  for (j=0;j<NITER;j++)
  {
    ia=(iswap=(index)) ^ c1[j];
    itmpl = ia & 0xffff;
    itmph = ia >> 16;
    ib=itmpl*itmpl+ ~(itmph*itmph);
    index=(lword) ^ (((ia = (ib >> 16) |
                      ((ib & 0xffff) << 16)) ^ c2[j])+itmpl*itmph);
    lword=iswap;
  }
  
  return ((double) index)*2.3283064365386963e-10;
}
#endif
#undef NITER
