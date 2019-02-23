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
isaac64.c: My random number generator for 64-bit machines.
By Bob Jenkins, 1996.  Public Domain.
------------------------------------------------------------------------------
*/
#include "config.h"
#include "isaac64.h"

#define mix(a, b, c, d, e, f, g, h)                                            \
  {                                                                            \
    a -= e;                                                                    \
    f ^= h >> 9;                                                               \
    h += a;                                                                    \
    b -= f;                                                                    \
    g ^= a << 9;                                                               \
    a += b;                                                                    \
    c -= g;                                                                    \
    h ^= b >> 23;                                                              \
    b += c;                                                                    \
    d -= h;                                                                    \
    a ^= c << 15;                                                              \
    c += d;                                                                    \
    e -= a;                                                                    \
    b ^= d >> 14;                                                              \
    d += e;                                                                    \
    f -= b;                                                                    \
    c ^= e << 20;                                                              \
    e += f;                                                                    \
    g -= c;                                                                    \
    d ^= f >> 17;                                                              \
    f += g;                                                                    \
    h -= d;                                                                    \
    e ^= g << 14;                                                              \
    g += h;                                                                    \
  }

void isaac64_init(struct isaac64_state *rng, ub4 seed) {
  ub8 *r, *m;
  ub8 a, b, c, d, e, f, g, h;
  ub4 i;

  rng->rngblocks = 0;

  rng->aa = (ub8)0;
  rng->bb = (ub8)0;
  rng->cc = (ub8)0;

  a = b = c = d = e = f = g = h = 0x9e3779b97f4a7c13LL; /* the golden ratio */

  r = rng->randrsl;
  m = rng->mm;

  for (i = 0; i < RANDSIZ; ++i)
    r[i] = (ub8)0;

  r[0] = seed;

  for (i = 0; i < 4; ++i) /* scramble it */
  {
    mix(a, b, c, d, e, f, g, h);
  }

  for (i = 0; i < RANDSIZ; i += 8) /* fill in m[] with messy stuff */
  {
    /* use all the information in the seed */
    a += r[i];
    b += r[i + 1];
    c += r[i + 2];
    d += r[i + 3];
    e += r[i + 4];
    f += r[i + 5];
    g += r[i + 6];
    h += r[i + 7];
    mix(a, b, c, d, e, f, g, h);
    m[i] = a;
    m[i + 1] = b;
    m[i + 2] = c;
    m[i + 3] = d;
    m[i + 4] = e;
    m[i + 5] = f;
    m[i + 6] = g;
    m[i + 7] = h;
  }

  /* do a second pass to make all of the seed affect all of m[] */
  for (i = 0; i < RANDSIZ; i += 8) {
    a += m[i];
    b += m[i + 1];
    c += m[i + 2];
    d += m[i + 3];
    e += m[i + 4];
    f += m[i + 5];
    g += m[i + 6];
    h += m[i + 7];
    mix(a, b, c, d, e, f, g, h);
    m[i] = a;
    m[i + 1] = b;
    m[i + 2] = c;
    m[i + 3] = d;
    m[i + 4] = e;
    m[i + 5] = f;
    m[i + 6] = g;
    m[i + 7] = h;
  }

  isaac64_generate(rng);  /* fill in the first set of results */
  rng->randcnt = RANDMAX; /* prepare to use the first set of results */
}
