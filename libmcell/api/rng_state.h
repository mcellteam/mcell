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

#ifndef API_RNG_STATE_H
#define API_RNG_STATE_H

#include "generated/gen_rng_state.h"
#include "api/api_common.h"

struct rng_state;

namespace MCell {
namespace API {

class RngState: public GenRngState {
public:
  RNG_STATE_CTOR()

  void check_semantics() const override {
    if (randslr.size() != RNG_SIZE) {
      throw ValueError(S("List ") + NAME_RANDSLR + " must have exactly " + std::to_string(RNG_SIZE) + " items.");
    }
    if (mm.size() != RNG_SIZE) {
      throw ValueError(S("List ") + NAME_MM + " must have exactly " + std::to_string(RNG_SIZE) + " items.");
    }
  }

  // internal, used for checkpointing
  RngState(const rng_state& rng);
};

} // namespace API
} // namespace MCell

#endif // API_RNG_STATE_H
