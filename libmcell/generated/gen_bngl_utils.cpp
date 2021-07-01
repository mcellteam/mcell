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
#include "gen_bngl_utils.h"

namespace MCell {
namespace API {

void define_pybinding_bngl_utils(py::module& m) {
  m.def_submodule("bngl_utils")
      .def("load_bngl_parameters", &bngl_utils::load_bngl_parameters, py::arg("file_name"), py::arg("parameter_overrides") = std::map<std::string, double>(), "Load parameters section from a BNGL file and return it as a dictionary name->value.\n- file_name: Path to the BNGL file to be loaded.\n\n- parameter_overrides: For each key k in the parameter_overrides, if it is defined in the BNGL's parameters section,\nits value is ignored and instead value parameter_overrides[k] is used.\n\n\n")
    ;
}

} // namespace API
} // namespace MCell

