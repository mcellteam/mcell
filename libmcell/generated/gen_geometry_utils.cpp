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
#include "gen_geometry_utils.h"
#include "api/geometry_object.h"

namespace MCell {
namespace API {

void define_pybinding_geometry_utils(py::module& m) {
  m.def_submodule("geometry_utils")
      .def("create_box", &geometry_utils::create_box, py::arg("name"), py::arg("edge_length"), "Creates a GeometryObject whose center is at (0, 0, 0).")
    ;
}

} // namespace API
} // namespace MCell

