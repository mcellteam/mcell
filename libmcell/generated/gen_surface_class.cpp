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
#include "gen_surface_class.h"
#include "../api/surface_class.h"
#include "../api/species.h"

namespace MCell {
namespace API {

void GenSurfaceClass::check_semantics() const {
  if (!is_set(name)) {
    throw ValueError("Parameter 'name' must be set.");
  }
}

bool GenSurfaceClass::__eq__(const GenSurfaceClass& other) const {
  return
    name == other.name &&
    name == other.name &&
    reflective->__eq__(*other.reflective) &&
    transparent->__eq__(*other.transparent) &&
    absorptive->__eq__(*other.absorptive);
}

void GenSurfaceClass::set_initialized() {
  if (is_set(reflective)) {
    reflective->set_initialized();
  }
  if (is_set(transparent)) {
    transparent->set_initialized();
  }
  if (is_set(absorptive)) {
    absorptive->set_initialized();
  }
  initialized = true;
}

std::string GenSurfaceClass::to_str(const std::string ind) const {
  std::stringstream ss;
  ss << get_object_name() << ": " <<
      "name=" << name << ", " <<
      "\n" << ind + "  " << "reflective=" << "(" << ((reflective != nullptr) ? reflective->to_str(ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "transparent=" << "(" << ((transparent != nullptr) ? transparent->to_str(ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "absorptive=" << "(" << ((absorptive != nullptr) ? absorptive->to_str(ind + "  ") : "null" ) << ")";
  return ss.str();
}

py::class_<SurfaceClass> define_pybinding_SurfaceClass(py::module& m) {
  return py::class_<SurfaceClass, std::shared_ptr<SurfaceClass>>(m, "SurfaceClass")
      .def(
          py::init<
            const std::string&,
            std::shared_ptr<Species>,
            std::shared_ptr<Species>,
            std::shared_ptr<Species>
          >(),
          py::arg("name"),
          py::arg("reflective") = nullptr,
          py::arg("transparent") = nullptr,
          py::arg("absorptive") = nullptr
      )
      .def("check_semantics", &SurfaceClass::check_semantics)
      .def("__str__", &SurfaceClass::to_str, py::arg("ind") = std::string(""))
      .def("dump", &SurfaceClass::dump)
      .def_property("name", &SurfaceClass::get_name, &SurfaceClass::set_name)
      .def_property("reflective", &SurfaceClass::get_reflective, &SurfaceClass::set_reflective)
      .def_property("transparent", &SurfaceClass::get_transparent, &SurfaceClass::set_transparent)
      .def_property("absorptive", &SurfaceClass::get_absorptive, &SurfaceClass::set_absorptive)
    ;
}

} // namespace API
} // namespace MCell

