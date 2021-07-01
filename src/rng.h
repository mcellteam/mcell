/******************************************************************************
 *
 * Copyright (C) 2006-2017 by
 * The Salk Institute for Biological Studies and
 * Pittsburgh Supercomputing Center, Carnegie Mellon University
 *
 * Use of this source code is governed by an MIT-style
 * license that can be found in the LICENSE file or at
 * https://opensource.org/licenses/MIT.
 *
******************************************************************************/

#pragma once

#include "debug_config.h"

#define ONE_OVER_2_TO_THE_33RD 1.16415321826934814453125e-10

#if defined(USE_MINIMAL_RNG)
#include "minrng.h"
#define rng_state mrng_state

#define rng_init(x, y) mrng_init((x), (y))
#define rng_dbl(x) mrng_dbl32((x))
#define rng_uint(x) mrng_uint32((x))

#else
/*******************ISAAC64*********************/
#include "isaac64.h"

#define rng_state isaac64_state

#define rng_uses(x)                                                            \
  ((RANDMAX *((x)->rngblocks - 1)) + (long long)(RANDMAX - (x)->randcnt))
#define rng_init(x, y) isaac64_init((x), (y))

#if !defined(NDEBUG) || defined(DEBUG_RNG_CALLS)
// we need functions to be able to dump the random number gen. info
static double rng_dbl(struct rng_state *rng);
static unsigned int rng_uint(struct rng_state *rng);
#else
#define rng_dbl(x) isaac64_dbl32((x))
#define rng_uint(x) isaac64_uint32((x))
#endif
/***********************************************/

#endif

#define rng_open_dbl(x) (rng_dbl(x) + ONE_OVER_2_TO_THE_33RD)

static double rng_gauss(struct rng_state *rng);

#include "rng.c"

