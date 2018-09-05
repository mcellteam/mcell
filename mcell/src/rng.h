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

#pragma once

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
#define rng_dbl(x) isaac64_dbl32((x))
#define rng_uint(x) isaac64_uint32((x))
/***********************************************/

#endif

#define rng_open_dbl(x) (rng_dbl(x) + ONE_OVER_2_TO_THE_33RD)

double rng_gauss(struct rng_state *rng);
