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
#include "gen_surface_region.h"
#include "../api/surface_region.h"
#include "../api/geometry_object.h"
#include "../api/region.h"

namespace MCell {
namespace API {

void GenSurfaceRegion::check_semantics() const {
  if (!is_set(name)) {
    throw ValueError("Parameter 'name' must be set.");
  }
  if (!is_set(element_connections)) {
    throw ValueError("Parameter 'element_connections' must be set.");
  }
}

bool GenSurfaceRegion::__eq__(const GenSurfaceRegion& other) const {
  return
    name == other.name &&
    name == other.name &&
    element_connections == other.element_connections &&
    parent->__eq__(*other.parent) &&
    node_type == other.node_type &&
    left_node->__eq__(*other.left_node) &&
    right_node->__eq__(*other.right_node);
}

void GenSurfaceRegion::set_initialized() {
  if (is_set(parent)) {
    parent->set_initialized();
  }
  if (is_set(left_node)) {
    left_node->set_initialized();
  }
  if (is_set(right_node)) {
    right_node->set_initialized();
  }
  initialized = true;
}

std::string GenSurfaceRegion::to_str(const std::string ind) const {
  std::stringstream ss;
  ss << get_object_name() << ": " <<
      "name=" << name << ", " <<
      "element_connections=" << vec_nonptr_to_str(element_connections, ind + "  ") << ", " <<
      "\n" << ind + "  " << "parent=" << "(" << ((parent != nullptr) ? parent->to_str(ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "node_type=" << node_type << ", " <<
      "\n" << ind + "  " << "left_node=" << "(" << ((left_node != nullptr) ? left_node->to_str(ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "right_node=" << "(" << ((right_node != nullptr) ? right_node->to_str(ind + "  ") : "null" ) << ")";
  return ss.str();
}

py::class_<SurfaceRegion> define_pybinding_SurfaceRegion(py::module& m) {
  return py::class_<SurfaceRegion, Region, std::shared_ptr<SurfaceRegion>>(m, "SurfaceRegion")
      .def(
          py::init<
            const std::string&,
            const std::vector<int>,
            std::shared_ptr<GeometryObject>,
            const RegionNodeType,
            std::shared_ptr<Region>,
            std::shared_ptr<Region>
          >(),
          py::arg("name"),
          py::arg("element_connections"),
          py::arg("parent") = nullptr,
          py::arg("node_type") = RegionNodeType::Unset,
          py::arg("left_node") = nullptr,
          py::arg("right_node") = nullptr
      )
      .def("check_semantics", &SurfaceRegion::check_semantics)
      .def("__str__", &SurfaceRegion::to_str, py::arg("ind") = std::string(""))
      .def("dump", &SurfaceRegion::dump)
      .def_property("name", &SurfaceRegion::get_name, &SurfaceRegion::set_name)
      .def_property("element_connections", &SurfaceRegion::get_element_connections, &SurfaceRegion::set_element_connections)
      .def_property("parent", &SurfaceRegion::get_parent, &SurfaceRegion::set_parent)
    ;
}

} // namespace API
} // namespace MCell

