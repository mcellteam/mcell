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
#include "gen_introspection.h"
#include "api/introspection.h"
#include "api/complex.h"
#include "api/geometry_object.h"
#include "api/molecule.h"
#include "api/wall.h"

namespace MCell {
namespace API {

bool GenIntrospection::__eq__(const Introspection& other) const {
  return
true ;
}

bool GenIntrospection::eq_nonarray_attributes(const Introspection& other, const bool ignore_name) const {
  return
true ;
}

std::string GenIntrospection::to_str(const std::string ind) const {
  std::stringstream ss;
  ss << "Introspection";
  return ss.str();
}

py::class_<Introspection> define_pybinding_Introspection(py::module& m) {
  return py::class_<Introspection, std::shared_ptr<Introspection>>(m, "Introspection", "This class is used only as a base class to Model, it is not provided through API. Provides methods to introspect simulation state.")
      .def(
          py::init<
          >()
      )
      .def("__str__", &Introspection::to_str, py::arg("ind") = std::string(""))
      .def("__eq__", &Introspection::__eq__, py::arg("other"))
      .def("get_molecule_ids", &Introspection::get_molecule_ids, py::arg("pattern") = nullptr, "Returns a list of ids of molecules.\nIf the arguments pattern is not set, the list of all molecule ids is returned.  \nIf the argument pattern is set, the list of all molecule ids whose species match \nthe pattern is returned. Matching of patterns with compartments works exactly in the \nsame was as in observables.\n")
      .def("get_molecule", &Introspection::get_molecule, py::arg("id"), "Returns a molecule from the simulated environment, None if the molecule does not exist")
      .def("get_species_name", &Introspection::get_species_name, py::arg("species_id"), "Returns a string representing canonical species name in the BNGL format. \n")
      .def("get_vertex", &Introspection::get_vertex, py::arg("object"), py::arg("vertex_index"), "Returns coordinates of a vertex.")
      .def("get_wall", &Introspection::get_wall, py::arg("object"), py::arg("wall_index"), "Returns information about a wall belonging to a given object.")
      .def("get_vertex_unit_normal", &Introspection::get_vertex_unit_normal, py::arg("object"), py::arg("vertex_index"), "Returns sum of all wall normals that use this vertex converted to a unit vector of length 1um.\nThis represents the unit vector pointing outwards from the vertex.\n")
      .def("get_wall_unit_normal", &Introspection::get_wall_unit_normal, py::arg("object"), py::arg("wall_index"), "Returns wall normal converted to a unit vector of length 1um.")
      .def("dump", &Introspection::dump)
    ;
}

} // namespace API
} // namespace MCell

