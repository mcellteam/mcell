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

#include <inttypes.h>

typedef uint32_t ub4;

#define rot(x, k) (((x) << (k)) | ((x) >> (32 - (k))))

struct mrng_state {
  ub4 a;
  ub4 b;
  ub4 c;
  ub4 d;
};

ub4 mrng_generate(struct mrng_state *x);
void mrng_init(struct mrng_state *x, ub4 seed);
#define mrng_uint32(rng) (mrng_generate(rng))

#define mrng_dbl32(rng) (DBL32 *(double)mrng_uint32(rng))
