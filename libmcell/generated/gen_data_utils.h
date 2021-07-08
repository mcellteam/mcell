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

#ifndef API_GEN_DATA_UTILS_H
#define API_GEN_DATA_UTILS_H

#include "api/api_common.h"

namespace MCell {
namespace API {

class PythonExportContext;

namespace data_utils {

std::vector<std::vector<double>> load_dat_file(const std::string& file_name);

} // namespace data_utils

void define_pybinding_data_utils(py::module& m);
} // namespace API
} // namespace MCell

#endif // API_GEN_DATA_UTILS_H
