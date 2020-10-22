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
#include "libs/pybind11/include/pybind11/stl.h"
#include "gen_surface_class.h"
#include "../api/surface_class.h"
#include "../api/species.h"
#include "../api/surface_property.h"

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
    vec_ptr_eq(properties, other.properties) &&
    type == other.type &&
    (
      (affected_species != nullptr) ?
        ( (other.affected_species != nullptr) ?
          (affected_species->__eq__(*other.affected_species)) : 
          false
        ) :
        ( (other.affected_species != nullptr) ?
          false :
          true
        )
     )  &&
    orientation == other.orientation;
}

void GenSurfaceClass::set_initialized() {
  vec_set_initialized(properties);
  if (is_set(affected_species)) {
    affected_species->set_initialized();
  }
  initialized = true;
}

void GenSurfaceClass::set_all_attributes_as_default_or_unset() {
  class_name = "SurfaceClass";
  name = STR_UNSET;
  properties = std::vector<std::shared_ptr<SurfaceProperty>>();
  type = SurfacePropertyType::UNSET;
  affected_species = nullptr;
  orientation = Orientation::NOT_SET;
}

std::string GenSurfaceClass::to_str(const std::string ind) const {
  std::stringstream ss;
  ss << get_object_name() << ": " <<
      "name=" << name << ", " <<
      "\n" << ind + "  " << "properties=" << vec_ptr_to_str(properties, ind + "  ") << ", " << "\n" << ind + "  " <<
      "type=" << type << ", " <<
      "\n" << ind + "  " << "affected_species=" << "(" << ((affected_species != nullptr) ? affected_species->to_str(ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "orientation=" << orientation;
  return ss.str();
}

py::class_<SurfaceClass> define_pybinding_SurfaceClass(py::module& m) {
  return py::class_<SurfaceClass, SurfaceProperty, std::shared_ptr<SurfaceClass>>(m, "SurfaceClass")
      .def(
          py::init<
            const std::string&,
            const std::vector<std::shared_ptr<SurfaceProperty>>,
            const SurfacePropertyType,
            std::shared_ptr<Species>,
            const Orientation
          >(),
          py::arg("name"),
          py::arg("properties") = std::vector<std::shared_ptr<SurfaceProperty>>(),
          py::arg("type") = SurfacePropertyType::UNSET,
          py::arg("affected_species") = nullptr,
          py::arg("orientation") = Orientation::NOT_SET
      )
      .def("check_semantics", &SurfaceClass::check_semantics)
      .def("__str__", &SurfaceClass::to_str, py::arg("ind") = std::string(""))
      .def("__repr__", &SurfaceClass::to_str, py::arg("ind") = std::string(""))
      .def("dump", &SurfaceClass::dump)
      .def_property("name", &SurfaceClass::get_name, &SurfaceClass::set_name)
      .def_property("properties", &SurfaceClass::get_properties, &SurfaceClass::set_properties)
    ;
}

} // namespace API
} // namespace MCell

