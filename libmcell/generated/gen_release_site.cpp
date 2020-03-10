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

// This file was generated automatically on 03/10/2020, 10:52 from 'data_classes.yaml'

#include <sstream>
#include "../api/mcell.h"
namespace MCell {
namespace API {

SemRes GenReleaseSite::check_semantics(std::ostream& out) const{
  if (!is_set(name)) {
    out << get_object_name() << ": Parameter 'name' must be set.\n";
    return SemRes::ERROR;
  }
  if (!is_set(shape)) {
    out << get_object_name() << ": Parameter 'shape' must be set.\n";
    return SemRes::ERROR;
  }
  if (!is_set(molecule)) {
    out << get_object_name() << ": Parameter 'molecule' must be set.\n";
    return SemRes::ERROR;
  }
  return SemRes::OK;
}

std::string GenReleaseSite::to_str() const{
  std::stringstream ss;
  ss << get_object_name() << ": " <<
      "name=" << name << ", " <<
      "shape=" << shape << ", " <<
      "molecule=" << ((molecule != nullptr) ? molecule->to_str() : "null" ) << ", " <<
      "location=" << location << ", " <<
      "site_diameter=" << site_diameter << ", " <<
      "site_radius=" << site_radius << ", " <<
      "release_probability=" << release_probability;
  return ss.str();
}

void define_binding_ReleaseSite(py::module& m) {
  py::class_<ReleaseSite>(m, "ReleaseSite")
      .def(
          py::init<
            const std::string&,
            const std::string&,
            const Species*,
            const Vec3&,
            const float_t,
            const float_t,
            const float_t
          >(),
          py::arg("name"),
          py::arg("shape"),
          py::arg("molecule"),
          py::arg("location") = VEC3_UNSET,
          py::arg("site_diameter") = FLT_UNSET,
          py::arg("site_radius") = FLT_UNSET,
          py::arg("release_probability") = FLT_UNSET
        )
      .def("check_semantics", &ReleaseSite::check_semantics_cerr)
      .def("__str__", &ReleaseSite::to_str)
      .def("dump", &ReleaseSite::dump)
      .def_property("name", &ReleaseSite::get_name, &ReleaseSite::set_name)
      .def_property("shape", &ReleaseSite::get_name, &ReleaseSite::set_name)
      .def_property("molecule", &ReleaseSite::get_name, &ReleaseSite::set_name)
      .def_property("location", &ReleaseSite::get_name, &ReleaseSite::set_name)
      .def_property("site_diameter", &ReleaseSite::get_name, &ReleaseSite::set_name)
      .def_property("site_radius", &ReleaseSite::get_name, &ReleaseSite::set_name)
      .def_property("release_probability", &ReleaseSite::get_name, &ReleaseSite::set_name)
    ;
}

} // namespace API
} // namespace MCell

