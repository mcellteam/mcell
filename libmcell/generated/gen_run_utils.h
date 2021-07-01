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

#ifndef API_GEN_RUN_UTILS_H
#define API_GEN_RUN_UTILS_H

#include "api/api_common.h"

namespace MCell {
namespace API {

class PythonExportContext;

namespace run_utils {

std::string get_last_checkpoint_dir(const int seed);
std::vector<std::string> remove_cwd(const std::vector<std::string> paths);

} // namespace run_utils

void define_pybinding_run_utils(py::module& m);
} // namespace API
} // namespace MCell

#endif // API_GEN_RUN_UTILS_H
