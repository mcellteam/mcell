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

/*
------------------------------------------------------------------------------
isaac64.h: definitions for a random number generator
Bob Jenkins, 1996, Public Domain

Modified for modularity by Tom Bartol and Rex Kerr
------------------------------------------------------------------------------
*/

#pragma once

#include <inttypes.h>

#define RANDSIZL (8)
#define RANDSIZ (1 << RANDSIZL)
#define RANDMAX (2 * RANDSIZ)

typedef unsigned long long ub8;
typedef uint32_t ub4;
typedef unsigned short int ub2;
typedef unsigned char ub1;

#define DBL32 (2.3283064365386962890625e-10)
#define DBL53 (1.1102230246251565404236316680908203125e-16)
#define DBL64 (5.42101086242752217003726400434970855712890625e-20)
#define MSK53 0x001FFFFFFFFFFFFFLL

struct isaac64_state {
  unsigned int randcnt;
  ub8 aa;
  ub8 bb;
  ub8 cc;
  ub8 randrsl[RANDSIZ];
  ub8 mm[RANDSIZ];
  ub8 rngblocks;
};

static void isaac64_init(struct isaac64_state *rng, ub4 seed);

static void isaac64_generate(struct isaac64_state *rng);

/*
------------------------------------------------------------------------------
Macros to get individual random numbers
------------------------------------------------------------------------------
*/

#define isaac64_uint32(rng)                                                    \
  ((rng)->randcnt > 0 ? (*(((ub4 *)((rng)->randrsl)) + ((rng)->randcnt -= 1)))       \
                    : (isaac64_generate(rng), (rng)->randcnt = RANDMAX - 1,      \
                       *(((ub4 *)((rng)->randrsl)) + (rng)->randcnt)))

#define isaac64_uint64(rng)                                                    \
  ((rng)->randcnt > 1                                                            \
       ? (*((ub8 *)(((ub4 *)((rng)->randrsl)) + ((rng)->randcnt -= 2))))           \
       : (isaac64_generate((rng)), (rng)->randcnt = RANDMAX - 2,                   \
          *((ub8 *)(((ub4 *)((rng)->randrsl)) + (rng)->randcnt))))

#define isaac64_dbl32(rng)                                                     \
  ((rng)->randcnt > 0                                                            \
       ? (DBL32 *(*(((ub4 *)((rng)->randrsl)) + ((rng)->randcnt -= 1))))           \
       : (isaac64_generate(rng), (rng)->randcnt = RANDMAX - 1,                   \
          DBL32 * (*(((ub4 *)((rng)->randrsl)) + (rng)->randcnt))))

#define isaac64_dbl53(rng)                                                     \
  ((rng)->randcnt > 1                                                            \
       ? (DBL53 *(                                                             \
             (*((ub8 *)(((ub4 *)((rng)->randrsl)) + ((rng)->randcnt -= 2)))) >>    \
             11))                                                              \
       : (isaac64_generate(rng), (rng)->randcnt = RANDMAX - 2,                   \
          DBL53 *                                                              \
              ((*((ub8 *)(((ub4 *)((rng)->randrsl)) + (rng)->randcnt))) >> 11)))

#define isaac64_dbl64(rng)                                                     \
  ((rng)->randcnt > 1                                                            \
       ? (DBL64 *(*((ub8 *)(((ub4 *)((rng)->randrsl)) + ((rng)->randcnt -= 2)))))  \
       : (isaac64_generate(rng), (rng)->randcnt = RANDMAX - 2,                   \
          DBL64 * (*((ub8 *)(((ub4 *)((rng)->randrsl)) + (rng)->randcnt)))))

#include "isaac64.c"
