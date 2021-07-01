/******************************************************************************
 *
 * Copyright (C) 2020 by
 * The Salk Institute for Biological Studies
 *
 * Use of this source code is governed by an MIT-style
 * license that can be found in the LICENSE file or at
 * https://opensource.org/licenses/MIT.
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
