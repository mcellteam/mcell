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

#include "generated/gen_geometry_utils.h"

#include "bng/bng.h"

#include "api/geometry_object.h"

using namespace std;

namespace MCell {
namespace API {

namespace bngl_utils {

std::map<std::string, double> load_bngl_parameters(
    const std::string& file_name,
    const std::map<std::string, double>& parameter_overrides
) {

  BNG::BNGData bng_data;

  int num_errors = BNG::parse_bngl_file(file_name, bng_data, parameter_overrides);
  if (num_errors != 0) {
    throw RuntimeError("Could not parse BNGL file " + file_name + ".");
  }

  return bng_data.get_parameters();
}

} // namespace bngl_utils

} // namespace API
} // namespace MCell
