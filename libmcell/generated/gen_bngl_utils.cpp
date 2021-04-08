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

#include <sstream>
#include "api/pybind11_stl_include.h"
#include "api/python_export_utils.h"
#include "gen_bngl_utils.h"

namespace MCell {
namespace API {

void define_pybinding_bngl_utils(py::module& m) {
  m.def_submodule("bngl_utils")
      .def("load_bngl_parameters", &bngl_utils::load_bngl_parameters, py::arg("file_name"), py::arg("parameter_overrides") = std::map<std::string, float_t>(), "Load parameters section from a BNGL file and return it as a dictionary name->value.\n- file_name: Path to the BNGL file to be loaded.\n\n- parameter_overrides: For each key k in the parameter_overrides, if it is defined in the BNGL's parameters section,\nits value is ignored and instead value parameter_overrides[k] is used.\n\n\n")
    ;
}

} // namespace API
} // namespace MCell

