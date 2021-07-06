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

#include <sstream>
#include "api/pybind11_stl_include.h"
#include "api/python_export_utils.h"
#include "gen_data_utils.h"

namespace MCell {
namespace API {

void define_pybinding_data_utils(py::module& m) {
  m.def_submodule("data_utils")
      .def("load_dat_file", &data_utils::load_dat_file, py::arg("file_name"), "Loads a two-column file where the first column is usually time and the second is a \nfloating point value. Returns a two-column list. \nCan be used to load a file with variable rate constants. \n\n- file_name: Path to the .dat file to be loaded.\n\n")
    ;
}

} // namespace API
} // namespace MCell

