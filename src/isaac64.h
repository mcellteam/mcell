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

void isaac64_init(struct isaac64_state *rng, ub4 seed);

#define ind(mm, x) (*(ub8 *)((ub1 *)(mm) + ((x) & ((RANDSIZ - 1) << 3))))

#define rngstep(mix, a, b, mm, m, m2, r, x)                                    \
  {                                                                            \
    x = *m;                                                                    \
    a = (mix) + *(m2++);                                                       \
    *(m++) = y = ind(mm, x) + a + b;                                           \
    *(r++) = b = ind(mm, y >> RANDSIZL) + x;                                   \
  }

static inline void isaac64_generate(struct isaac64_state *rng) {
  register ub8 a, b, x, y, *m, *m2, *r, *mend;

  m = rng->mm;
  r = rng->randrsl;
  a = rng->aa;
  b = rng->bb + (++rng->cc);
  for (m = rng->mm, mend = m2 = m + (RANDSIZ / 2); m < mend;) {
    rngstep(~(a ^ (a << 21)), a, b, rng->mm, m, m2, r, x);
    rngstep(a ^ (a >> 5), a, b, rng->mm, m, m2, r, x);
    rngstep(a ^ (a << 12), a, b, rng->mm, m, m2, r, x);
    rngstep(a ^ (a >> 33), a, b, rng->mm, m, m2, r, x);
  }
  for (m2 = rng->mm; m2 < mend;) {
    rngstep(~(a ^ (a << 21)), a, b, rng->mm, m, m2, r, x);
    rngstep(a ^ (a >> 5), a, b, rng->mm, m, m2, r, x);
    rngstep(a ^ (a << 12), a, b, rng->mm, m, m2, r, x);
    rngstep(a ^ (a >> 33), a, b, rng->mm, m, m2, r, x);
  }
  rng->bb = b;
  rng->aa = a;
  ++rng->rngblocks;
}
/*
------------------------------------------------------------------------------
Macros to get individual random numbers
------------------------------------------------------------------------------
*/

#define isaac64_uint32(rng)                                                    \
  (rng->randcnt > 0 ? (*(((ub4 *)(rng->randrsl)) + (rng->randcnt -= 1)))       \
                    : (isaac64_generate(rng), rng->randcnt = RANDMAX - 1,      \
                       *(((ub4 *)(rng->randrsl)) + rng->randcnt)))

#define isaac64_uint64(rng)                                                    \
  (rng->randcnt > 1                                                            \
       ? (*((ub8 *)(((ub4 *)(rng->randrsl)) + (rng->randcnt -= 2))))           \
       : (isaac64_generate(rng), rng->randcnt = RANDMAX - 2,                   \
          *((ub8 *)(((ub4 *)(rng->randrsl)) + rng->randcnt))))

#define isaac64_dbl32(rng)                                                     \
  (rng->randcnt > 0                                                            \
       ? (DBL32 *(*(((ub4 *)(rng->randrsl)) + (rng->randcnt -= 1))))           \
       : (isaac64_generate(rng), rng->randcnt = RANDMAX - 1,                   \
          DBL32 * (*(((ub4 *)(rng->randrsl)) + rng->randcnt))))

#define isaac64_dbl53(rng)                                                     \
  (rng->randcnt > 1                                                            \
       ? (DBL53 *(                                                             \
             (*((ub8 *)(((ub4 *)(rng->randrsl)) + (rng->randcnt -= 2)))) >>    \
             11))                                                              \
       : (isaac64_generate(rng), rng->randcnt = RANDMAX - 2,                   \
          DBL53 *                                                              \
              ((*((ub8 *)(((ub4 *)(rng->randrsl)) + rng->randcnt))) >> 11)))

#define isaac64_dbl64(rng)                                                     \
  (rng->randcnt > 1                                                            \
       ? (DBL64 *(*((ub8 *)(((ub4 *)(rng->randrsl)) + (rng->randcnt -= 2)))))  \
       : (isaac64_generate(rng), rng->randcnt = RANDMAX - 2,                   \
          DBL64 * (*((ub8 *)(((ub4 *)(rng->randrsl)) + rng->randcnt)))))
