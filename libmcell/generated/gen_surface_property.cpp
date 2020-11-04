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
#include "gen_surface_property.h"
#include "../api/surface_property.h"
#include "../api/complex.h"

namespace MCell {
namespace API {

void GenSurfaceProperty::check_semantics() const {
}

bool GenSurfaceProperty::__eq__(const GenSurfaceProperty& other) const {
  return
    name == other.name &&
    type == other.type &&
    (
      (affected_complex_pattern != nullptr) ?
        ( (other.affected_complex_pattern != nullptr) ?
          (affected_complex_pattern->__eq__(*other.affected_complex_pattern)) : 
          false
        ) :
        ( (other.affected_complex_pattern != nullptr) ?
          false :
          true
        )
     )  &&
    clamp_concentration == other.clamp_concentration;
}

void GenSurfaceProperty::set_initialized() {
  if (is_set(affected_complex_pattern)) {
    affected_complex_pattern->set_initialized();
  }
  initialized = true;
}

void GenSurfaceProperty::set_all_attributes_as_default_or_unset() {
  class_name = "SurfaceProperty";
  type = SurfacePropertyType::UNSET;
  affected_complex_pattern = nullptr;
  clamp_concentration = FLT_UNSET;
}

std::string GenSurfaceProperty::to_str(const std::string ind) const {
  std::stringstream ss;
  ss << get_object_name() << ": " <<
      "type=" << type << ", " <<
      "\n" << ind + "  " << "affected_complex_pattern=" << "(" << ((affected_complex_pattern != nullptr) ? affected_complex_pattern->to_str(ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "clamp_concentration=" << clamp_concentration;
  return ss.str();
}

py::class_<SurfaceProperty> define_pybinding_SurfaceProperty(py::module& m) {
  return py::class_<SurfaceProperty, std::shared_ptr<SurfaceProperty>>(m, "SurfaceProperty")
      .def(
          py::init<
            const SurfacePropertyType,
            std::shared_ptr<Complex>,
            const float_t
          >(),
          py::arg("type") = SurfacePropertyType::UNSET,
          py::arg("affected_complex_pattern") = nullptr,
          py::arg("clamp_concentration") = FLT_UNSET
      )
      .def("check_semantics", &SurfaceProperty::check_semantics)
      .def("__str__", &SurfaceProperty::to_str, py::arg("ind") = std::string(""))
      .def("__repr__", &SurfaceProperty::to_str, py::arg("ind") = std::string(""))
      .def("dump", &SurfaceProperty::dump)
      .def_property("type", &SurfaceProperty::get_type, &SurfaceProperty::set_type)
      .def_property("affected_complex_pattern", &SurfaceProperty::get_affected_complex_pattern, &SurfaceProperty::set_affected_complex_pattern)
      .def_property("clamp_concentration", &SurfaceProperty::get_clamp_concentration, &SurfaceProperty::set_clamp_concentration)
    ;
}

} // namespace API
} // namespace MCell

