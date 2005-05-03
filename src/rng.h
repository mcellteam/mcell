#ifndef RNG_H
#define RNG_H

/*#define USE_RAN4*/

#ifdef USE_RAN4

/*********************RAN4**********************/
#include "ran4.h"

#define rng_state ran4_state

#define rng_init(x,y) ran4_init((x),(y))
#define rng_dbl(x) ran4_dbl32((x))
#define rng_uint(x) ran4_uint32((x))
/***********************************************/

#else

/*******************ISAAC64*********************/
#include "isaac64.h"

#define rng_state isaac64_state

#define rng_init(x,y) isaac64_init((x),(y))
#define rng_dbl(x) isaac64_dbl32((x))
#define rng_uint(x) isaac64_uint32((x))
/***********************************************/

#endif

double rng_gauss(struct rng_state *rng);

#endif

