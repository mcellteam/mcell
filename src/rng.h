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

#include <math.h>

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

//double rng_gauss(struct rng_state *rng);


#define R_VALUE (3.442619855899)
static const double SCALE_FACTOR = R_VALUE;
static const double RECIP_SCALE_FACTOR = 1.0 / R_VALUE;

/* Tabulated PDF at ends of strips */
extern const double YTAB[128];

/* Tabulated 'K' for quick out on strips (strip #0 at bottom, strips
 * 1...127 counting down from the top */
extern const unsigned long KTAB[128];

/* Tabulated 'W' - scale for output values which aren't in the tails.
 */
extern const double WTAB[128];

/*************************************************************************
rng_gauss:
  In:  struct rng_state *rng - uniform RNG state
  Out: Returns a Gaussian variate (mean 0, variance 1)
 *************************************************************************/
static inline double rng_gauss(struct rng_state *rng) {
  double x, y;
  double sign = 1.0;

  int npasses = 0;
  do {
    unsigned long bits = rng_uint(rng);
    unsigned long region, pos_within_region;
    ++npasses;

    /* Partition bits:
     *    - Bits 0...7: select a region under the curve
     *    - Bit 8:      sign bit
     *    - Bits 9...31 pick a point within the region
     * */
    sign = (bits & 0x80) ? -1.0 : 1.0;
    region = bits & 0x0000007f;
    pos_within_region = bits & 0xffffff00;

    /* Compute our X, and check if the X value lies entirely under the
     * curve for the chosen region */
    x = pos_within_region * WTAB[region];
    if (pos_within_region < KTAB[region])
      break;

    /* If we're in one of the 127 cheap regions */
    if (region != 0) {
      double yR, yB;
      yB = YTAB[region];
      yR = YTAB[region - 1] - yB;
      y = yB + yR * rng_dbl(rng);
    }

    /* If we're in the expensive region */
    else {
      x = SCALE_FACTOR - log1p(-rng_dbl(rng)) * RECIP_SCALE_FACTOR;
      y = exp(-SCALE_FACTOR * (x - 0.5 * SCALE_FACTOR)) * rng_dbl(rng);
    }
  } while (y >= exp(-0.5 * x * x));

  //printf("rng_gauss: %f\n", sign * x);
  return sign * x;
}

