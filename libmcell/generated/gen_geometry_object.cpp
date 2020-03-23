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

// This file was generated automatically on 03/23/2020, 15:47 from 'data_classes.yaml'

#include <sstream>
#include "../api/mcell.h"
namespace MCell {
namespace API {

SemRes GenGeometryObject::check_semantics(std::ostream& out) const{
  if (!is_set(name)) {
    out << get_object_name() << ": Parameter 'name' must be set.\n";
    return SemRes::ERROR;
  }
  return SemRes::OK;
}

std::string GenGeometryObject::to_str() const{
  std::stringstream ss;
  ss << get_object_name() << ": " <<
      "name=" << name;
  return ss.str();
}

void define_pybinding_GeometryObject(py::module& m) {
  py::class_<GeometryObject>(m, "GeometryObject")
      .def(
          py::init<
            const std::string&
          >(),
          py::arg("name")
        )
      .def("check_semantics", &GeometryObject::check_semantics_cerr)
      .def("__str__", &GeometryObject::to_str)
      .def("dump", &GeometryObject::dump)
      .def_property("name", &GeometryObject::get_name, &GeometryObject::set_name)
    ;
}

} // namespace API
} // namespace MCell

