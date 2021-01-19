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
