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
#include "gen_region.h"
#include "../api/region.h"
#include "../api/region.h"

namespace MCell {
namespace API {

void GenRegion::check_semantics() const {
}

bool GenRegion::__eq__(const GenRegion& other) const {
  return
    name == other.name &&
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

void GenRegion::set_initialized() {
  if (is_set(left_node)) {
    left_node->set_initialized();
  }
  if (is_set(right_node)) {
    right_node->set_initialized();
  }
  initialized = true;
}

std::string GenRegion::to_str(const std::string ind) const {
  std::stringstream ss;
  ss << get_object_name() << ": " <<
      "node_type=" << node_type << ", " <<
      "\n" << ind + "  " << "left_node=" << "(" << ((left_node != nullptr) ? left_node->to_str(ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "right_node=" << "(" << ((right_node != nullptr) ? right_node->to_str(ind + "  ") : "null" ) << ")";
  return ss.str();
}

py::class_<Region> define_pybinding_Region(py::module& m) {
  return py::class_<Region, std::shared_ptr<Region>>(m, "Region")
      .def(
          py::init<
            const RegionNodeType,
            std::shared_ptr<Region>,
            std::shared_ptr<Region>
          >(),
          py::arg("node_type") = RegionNodeType::UNSET,
          py::arg("left_node") = nullptr,
          py::arg("right_node") = nullptr
      )
      .def("check_semantics", &Region::check_semantics)
      .def("__str__", &Region::to_str, py::arg("ind") = std::string(""))
      .def("__add__", &Region::__add__, py::arg("other"))
      .def("__sub__", &Region::__sub__, py::arg("other"))
      .def("__mul__", &Region::__mul__, py::arg("other"))
      .def("dump", &Region::dump)
      .def_property("node_type", &Region::get_node_type, &Region::set_node_type)
      .def_property("left_node", &Region::get_left_node, &Region::set_left_node)
      .def_property("right_node", &Region::get_right_node, &Region::set_right_node)
    ;
}

} // namespace API
} // namespace MCell

