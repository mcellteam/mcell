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
#include "gen_volume_compartment.h"
#include "../api/volume_compartment.h"
#include "../api/geometry_object.h"
#include "../api/region.h"
#include "../api/volume_compartment.h"

namespace MCell {
namespace API {

void GenVolumeCompartment::check_semantics() const {
  if (!is_set(name)) {
    throw ValueError("Parameter 'name' must be set.");
  }
  if (!is_set(geometry_object)) {
    throw ValueError("Parameter 'geometry_object' must be set.");
  }
}

bool GenVolumeCompartment::__eq__(const GenVolumeCompartment& other) const {
  return
    name == other.name &&
    name == other.name &&
    (
      (geometry_object != nullptr) ?
        ( (other.geometry_object != nullptr) ?
          (geometry_object->__eq__(*other.geometry_object)) : 
          false
        ) :
        ( (other.geometry_object != nullptr) ?
          false :
          true
        )
     )  &&
    vec_ptr_eq(child_compartments, other.child_compartments) &&
    surface_compartment_name == other.surface_compartment_name;
}

void GenVolumeCompartment::set_initialized() {
  if (is_set(geometry_object)) {
    geometry_object->set_initialized();
  }
  vec_set_initialized(child_compartments);
  initialized = true;
}

void GenVolumeCompartment::set_all_attributes_as_default_or_unset() {
  class_name = "VolumeCompartment";
  name = STR_UNSET;
  geometry_object = nullptr;
  child_compartments = std::vector<std::shared_ptr<VolumeCompartment>>();
  surface_compartment_name = STR_UNSET;
}

std::string GenVolumeCompartment::to_str(const std::string ind) const {
  std::stringstream ss;
  ss << get_object_name() << ": " <<
      "name=" << name << ", " <<
      "\n" << ind + "  " << "geometry_object=" << "(" << ((geometry_object != nullptr) ? geometry_object->to_str(ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "child_compartments=" << vec_ptr_to_str(child_compartments, ind + "  ") << ", " << "\n" << ind + "  " <<
      "surface_compartment_name=" << surface_compartment_name;
  return ss.str();
}

py::class_<VolumeCompartment> define_pybinding_VolumeCompartment(py::module& m) {
  return py::class_<VolumeCompartment, std::shared_ptr<VolumeCompartment>>(m, "VolumeCompartment")
      .def(
          py::init<
            const std::string&,
            std::shared_ptr<GeometryObject>,
            const std::vector<std::shared_ptr<VolumeCompartment>>,
            const std::string&
          >(),
          py::arg("name"),
          py::arg("geometry_object"),
          py::arg("child_compartments") = std::vector<std::shared_ptr<VolumeCompartment>>(),
          py::arg("surface_compartment_name") = STR_UNSET
      )
      .def("check_semantics", &VolumeCompartment::check_semantics)
      .def("__str__", &VolumeCompartment::to_str, py::arg("ind") = std::string(""))
      .def("get_volume_compartment_region", &VolumeCompartment::get_volume_compartment_region)
      .def("dump", &VolumeCompartment::dump)
      .def_property("name", &VolumeCompartment::get_name, &VolumeCompartment::set_name)
      .def_property("geometry_object", &VolumeCompartment::get_geometry_object, &VolumeCompartment::set_geometry_object)
      .def_property("child_compartments", &VolumeCompartment::get_child_compartments, &VolumeCompartment::set_child_compartments)
      .def_property("surface_compartment_name", &VolumeCompartment::get_surface_compartment_name, &VolumeCompartment::set_surface_compartment_name)
    ;
}

} // namespace API
} // namespace MCell

