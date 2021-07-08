/******************************************************************************
 *
 * Copyright (C) 2021 by
 * The Salk Institute for Biological Studies
 *
 * Use of this source code is governed by an MIT-style
 * license that can be found in the LICENSE file or at
 * https://opensource.org/licenses/MIT.
 *
******************************************************************************/

#ifndef API_GEN_BNGL_UTILS_H
#define API_GEN_BNGL_UTILS_H

#include "api/api_common.h"

namespace MCell {
namespace API {

class PythonExportContext;

namespace bngl_utils {

std::map<std::string, double> load_bngl_parameters(const std::string& file_name, const std::map<std::string, double>& parameter_overrides = std::map<std::string, double>());

} // namespace bngl_utils

void define_pybinding_bngl_utils(py::module& m);
} // namespace API
} // namespace MCell

#endif // API_GEN_BNGL_UTILS_H
