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

#ifndef API_VIZ_OUTPUT_H
#define API_VIZ_OUTPUT_H

#include "generated/gen_viz_output.h"
#include "api/api_common.h"

namespace MCell {
namespace API {

class VizOutput: public GenVizOutput {
public:
  VIZ_OUTPUT_CTOR()

  void check_semantics() const override {
    GenVizOutput::check_semantics();

    if (every_n_timesteps < 0) {
      throw ValueError(
          S("The value of ") + NAME_EVERY_N_TIMESTEPS + " must not be less than 0.");
    }
  }
};

} // namespace API
} // namespace MCell

#endif // API_VIZ_OUTPUT_H
