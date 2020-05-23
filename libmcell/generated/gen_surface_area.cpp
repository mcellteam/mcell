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
#include <pybind11/stl.h>
#include "gen_surface_area.h"
#include "../api/surface_area.h"
#include "../api/geometry_object.h"
#include "../api/region.h"

namespace MCell {
namespace API {

void GenSurfaceArea::check_semantics() const {
  if (!is_set(name)) {
    throw ValueError("Parameter 'name' must be set.");
  }
  if (!is_set(element_connections)) {
    throw ValueError("Parameter 'element_connections' must be set.");
  }
}

bool GenSurfaceArea::__eq__(const GenSurfaceArea& other) const {
  return
    name == other.name &&
    name == other.name &&
    element_connections == other.element_connections &&
    parent->__eq__(*other.parent);
}

void GenSurfaceArea::set_initialized() {
  parent->set_initialized();
  initialized = true;
}

std::string GenSurfaceArea::to_str(const std::string ind) const {
  std::stringstream ss;
  ss << get_object_name() << ": " <<
      "name=" << name << ", " <<
      "element_connections=" << vec_nonptr_to_str(element_connections, ind + "  ") << ", " <<
      "\n" << ind + "  " << "parent=" << "(" << ((parent != nullptr) ? parent->to_str(ind + "  ") : "null" ) << ")";
  return ss.str();
}

py::class_<SurfaceArea> define_pybinding_SurfaceArea(py::module& m) {
  return py::class_<SurfaceArea, std::shared_ptr<SurfaceArea>>(m, "SurfaceArea")
      .def(
          py::init<
            const std::string&,
            const std::vector<int>,
            std::shared_ptr<GeometryObject>
          >(),
          py::arg("name"),
          py::arg("element_connections"),
          py::arg("parent") = nullptr
      )
      .def("check_semantics", &SurfaceArea::check_semantics)
      .def("__str__", &SurfaceArea::to_str, py::arg("ind") = std::string(""))
      .def("as_region", &SurfaceArea::as_region)
      .def("dump", &SurfaceArea::dump)
      .def_property("name", &SurfaceArea::get_name, &SurfaceArea::set_name)
      .def_property("element_connections", &SurfaceArea::get_element_connections, &SurfaceArea::set_element_connections)
      .def_property("parent", &SurfaceArea::get_parent, &SurfaceArea::set_parent)
    ;
}

} // namespace API
} // namespace MCell

