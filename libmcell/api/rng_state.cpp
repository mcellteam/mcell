/******************************************************************************
 *
 * Copyright (C) 2020 by
 * The Salk Institute for Biological Studies
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

#include "api/rng_state.h"

#include "rng.h"

using namespace std;

namespace MCell {
namespace API {

RngState::RngState(const rng_state& rng) {
  assert(RNG_SIZE == RANDSIZ);

  randcnt = rng.randcnt;
  aa = rng.aa;
  bb = rng.bb;
  cc = rng.cc;
  cc = rng.cc;
  std::copy(&rng.randrsl[0], &rng.randrsl[RANDSIZ], std::back_inserter(randslr));
  std::copy(&rng.mm[0], &rng.mm[RANDSIZ], std::back_inserter(mm));

  rngblocks = rng.rngblocks; // counter for simulation stats, does not affect simulation
}

} // namespace API
} // namespace MCell
