/******************************************************************************
 *
 * Copyright (C) 2020 by
 * The Salk Institute for Biological Studies and
 * Pittsburgh Supercomputing Center, Carnegie Mellon University
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
#include "gen_surface_region.h"
#include "../api/surface_region.h"

namespace MCell {
namespace API {

SemRes GenSurfaceRegion::check_semantics(std::ostream& out) const {
  if (!is_set(name)) {
    out << get_object_name() << ": Parameter 'name' must be set.\n";
    return SemRes::ERROR;
  }
  if (!is_set(element_connections)) {
    out << get_object_name() << ": Parameter 'element_connections' must be set.\n";
    return SemRes::ERROR;
  }
  return SemRes::OK;
}

std::string GenSurfaceRegion::to_str(const std::string ind) const {
  std::stringstream ss;
  ss << get_object_name() << ": " <<
      "name=" << name << ", " <<
      "element_connections=" << vec_nonptr_to_str(element_connections, ind + "  ");
  return ss.str();
}

py::class_<SurfaceRegion> define_pybinding_SurfaceRegion(py::module& m) {
  return py::class_<SurfaceRegion, std::shared_ptr<SurfaceRegion>>(m, "SurfaceRegion")
      .def(
          py::init<
            const std::string&,
            const std::vector<int>
          >(),
          py::arg("name"),
          py::arg("element_connections")
        )
      .def("check_semantics", &SurfaceRegion::check_semantics_cerr)
      .def("__str__", &SurfaceRegion::to_str, py::arg("ind") = std::string(""))
      .def("dump", &SurfaceRegion::dump)
      .def_property("name", &SurfaceRegion::get_name, &SurfaceRegion::set_name)
      .def_property("element_connections", &SurfaceRegion::get_element_connections, &SurfaceRegion::set_element_connections)
    ;
}

} // namespace API
} // namespace MCell

