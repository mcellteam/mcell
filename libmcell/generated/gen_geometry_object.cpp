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
#include "gen_geometry_object.h"
#include "../api/geometry_object.h"
#include "../api/region.h"
#include "../api/surface_class.h"
#include "../api/surface_region.h"

namespace MCell {
namespace API {

void GenGeometryObject::check_semantics() const {
  if (!is_set(name)) {
    throw ValueError("Parameter 'name' must be set.");
  }
  if (!is_set(vertex_list)) {
    throw ValueError("Parameter 'vertex_list' must be set.");
  }
  if (!is_set(element_connections)) {
    throw ValueError("Parameter 'element_connections' must be set.");
  }
}

bool GenGeometryObject::__eq__(const GenGeometryObject& other) const {
  return
    name == other.name &&
    name == other.name &&
    vertex_list == other.vertex_list &&
    element_connections == other.element_connections &&
    vec_ptr_eq(surface_regions, other.surface_regions) &&
    (
      (surface_class != nullptr) ?
        ( (other.surface_class != nullptr) ?
          (surface_class->__eq__(*other.surface_class)) : 
          false
        ) :
        ( (other.surface_class != nullptr) ?
          false :
          true
        )
     )  &&
    node_type == other.node_type &&
    (
      (left_node != nullptr) ?
        ( (other.left_node != nullptr) ?
          (left_node->__eq__(*other.left_node)) : 
          false
        ) :
        ( (other.left_node != nullptr) ?
          false :
          true
        )
     )  &&
    (
      (right_node != nullptr) ?
        ( (other.right_node != nullptr) ?
          (right_node->__eq__(*other.right_node)) : 
          false
        ) :
        ( (other.right_node != nullptr) ?
          false :
          true
        )
     ) ;
}

void GenGeometryObject::set_initialized() {
  vec_set_initialized(surface_regions);
  if (is_set(surface_class)) {
    surface_class->set_initialized();
  }
  if (is_set(left_node)) {
    left_node->set_initialized();
  }
  if (is_set(right_node)) {
    right_node->set_initialized();
  }
  initialized = true;
}

std::string GenGeometryObject::to_str(const std::string ind) const {
  std::stringstream ss;
  ss << get_object_name() << ": " <<
      "name=" << name << ", " <<
      "vertex_list=" << vec_nonptr_to_str(vertex_list, ind + "  ") << ", " <<
      "element_connections=" << vec_nonptr_to_str(element_connections, ind + "  ") << ", " <<
      "\n" << ind + "  " << "surface_regions=" << vec_ptr_to_str(surface_regions, ind + "  ") << ", " << "\n" << ind + "  " <<
      "surface_class=" << "(" << ((surface_class != nullptr) ? surface_class->to_str(ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "node_type=" << node_type << ", " <<
      "\n" << ind + "  " << "left_node=" << "(" << ((left_node != nullptr) ? left_node->to_str(ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "right_node=" << "(" << ((right_node != nullptr) ? right_node->to_str(ind + "  ") : "null" ) << ")";
  return ss.str();
}

py::class_<GeometryObject> define_pybinding_GeometryObject(py::module& m) {
  return py::class_<GeometryObject, Region, std::shared_ptr<GeometryObject>>(m, "GeometryObject")
      .def(
          py::init<
            const std::string&,
            const std::vector<std::vector<float_t>>,
            const std::vector<std::vector<int>>,
            const std::vector<std::shared_ptr<SurfaceRegion>>,
            std::shared_ptr<SurfaceClass>,
            const RegionNodeType,
            std::shared_ptr<Region>,
            std::shared_ptr<Region>
          >(),
          py::arg("name"),
          py::arg("vertex_list"),
          py::arg("element_connections"),
          py::arg("surface_regions") = std::vector<std::shared_ptr<SurfaceRegion>>(),
          py::arg("surface_class") = nullptr,
          py::arg("node_type") = RegionNodeType::UNSET,
          py::arg("left_node") = nullptr,
          py::arg("right_node") = nullptr
      )
      .def("check_semantics", &GeometryObject::check_semantics)
      .def("__str__", &GeometryObject::to_str, py::arg("ind") = std::string(""))
      .def("dump", &GeometryObject::dump)
      .def_property("name", &GeometryObject::get_name, &GeometryObject::set_name)
      .def_property("vertex_list", &GeometryObject::get_vertex_list, &GeometryObject::set_vertex_list)
      .def_property("element_connections", &GeometryObject::get_element_connections, &GeometryObject::set_element_connections)
      .def_property("surface_regions", &GeometryObject::get_surface_regions, &GeometryObject::set_surface_regions)
      .def_property("surface_class", &GeometryObject::get_surface_class, &GeometryObject::set_surface_class)
    ;
}

} // namespace API
} // namespace MCell

