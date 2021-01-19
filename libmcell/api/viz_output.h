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
