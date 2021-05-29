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
#include "api/pybind11_stl_include.h"
#include "api/python_export_utils.h"
#include "gen_region.h"
#include "api/region.h"
#include "api/region.h"

namespace MCell {
namespace API {

void GenRegion::check_semantics() const {
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

void GenRegion::set_all_attributes_as_default_or_unset() {
  class_name = "Region";
  node_type = RegionNodeType::UNSET;
  left_node = nullptr;
  right_node = nullptr;
}

std::shared_ptr<Region> GenRegion::copy_region() const {
  std::shared_ptr<Region> res = std::make_shared<Region>(DefaultCtorArgType());
  res->class_name = class_name;
  res->node_type = node_type;
  res->left_node = left_node;
  res->right_node = right_node;

  return res;
}

std::shared_ptr<Region> GenRegion::deepcopy_region(py::dict) const {
  std::shared_ptr<Region> res = std::make_shared<Region>(DefaultCtorArgType());
  res->class_name = class_name;
  res->node_type = node_type;
  res->left_node = is_set(left_node) ? left_node->deepcopy_region() : nullptr;
  res->right_node = is_set(right_node) ? right_node->deepcopy_region() : nullptr;

  return res;
}

bool GenRegion::__eq__(const Region& other) const {
  return
    node_type == other.node_type &&
    (
      (is_set(left_node)) ?
        (is_set(other.left_node) ?
          (left_node->__eq__(*other.left_node)) : 
          false
        ) :
        (is_set(other.left_node) ?
          false :
          true
        )
     )  &&
    (
      (is_set(right_node)) ?
        (is_set(other.right_node) ?
          (right_node->__eq__(*other.right_node)) : 
          false
        ) :
        (is_set(other.right_node) ?
          false :
          true
        )
     ) ;
}

bool GenRegion::eq_nonarray_attributes(const Region& other, const bool ignore_name) const {
  return
    node_type == other.node_type &&
    (
      (is_set(left_node)) ?
        (is_set(other.left_node) ?
          (left_node->__eq__(*other.left_node)) : 
          false
        ) :
        (is_set(other.left_node) ?
          false :
          true
        )
     )  &&
    (
      (is_set(right_node)) ?
        (is_set(other.right_node) ?
          (right_node->__eq__(*other.right_node)) : 
          false
        ) :
        (is_set(other.right_node) ?
          false :
          true
        )
     ) ;
}

std::string GenRegion::to_str(const bool all_details, const std::string ind) const {
  std::stringstream ss;
  ss << get_object_name() << ": " <<
      "node_type=" << node_type << ", " <<
      "\n" << ind + "  " << "left_node=" << "(" << ((left_node != nullptr) ? left_node->to_str(all_details, ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "right_node=" << "(" << ((right_node != nullptr) ? right_node->to_str(all_details, ind + "  ") : "null" ) << ")";
  return ss.str();
}

py::class_<Region> define_pybinding_Region(py::module& m) {
  return py::class_<Region, std::shared_ptr<Region>>(m, "Region", "Represents region construted from 1 or more multiple, usually unnamed?")
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
      .def("__copy__", &Region::copy_region)
      .def("__deepcopy__", &Region::deepcopy_region, py::arg("memo"))
      .def("__str__", &Region::to_str, py::arg("all_details") = false, py::arg("ind") = std::string(""))
      .def("__eq__", &Region::__eq__, py::arg("other"))
      .def("__add__", &Region::__add__, py::arg("other"), "Computes union of two regions, use with Python operator '+'.\n- other\n")
      .def("__sub__", &Region::__sub__, py::arg("other"), "Computes difference of two regions, use with Python operator '-'.\n- other\n")
      .def("__mul__", &Region::__mul__, py::arg("other"), "Computes intersection of two regions, use with Python operator '*'.\n- other\n")
      .def("dump", &Region::dump)
      .def_property("node_type", &Region::get_node_type, &Region::set_node_type, "When this values is LeafGeometryObject, then this object is of class GeometryObject,\nwhen LeafSurfaceRegion, then it is of class SurfaceRegion.\n")
      .def_property("left_node", &Region::get_left_node, &Region::set_left_node, "Internal, do not use. When node_type is not Leaf, this is the left operand")
      .def_property("right_node", &Region::get_right_node, &Region::set_right_node, "Internal, do not use. When node_type is not Leaf, this is the right operand")
    ;
}

std::string GenRegion::export_to_python(std::ostream& out, PythonExportContext& ctx) {
  if (!export_even_if_already_exported() && ctx.already_exported(this)) {
    return ctx.get_exported_name(this);
  }
  std::string exported_name = "region_" + std::to_string(ctx.postinc_counter("region"));
  if (!export_even_if_already_exported()) {
    ctx.add_exported(this, exported_name);
  }

  bool str_export = export_as_string_without_newlines();
  std::string nl = "";
  std::string ind = " ";
  std::stringstream ss;
  if (!str_export) {
    nl = "\n";
    ind = "    ";
    ss << exported_name << " = ";
  }
  ss << "m.Region(" << nl;
  if (node_type != RegionNodeType::UNSET) {
    ss << ind << "node_type = " << node_type << "," << nl;
  }
  if (is_set(left_node)) {
    ss << ind << "left_node = " << left_node->export_to_python(out, ctx) << "," << nl;
  }
  if (is_set(right_node)) {
    ss << ind << "right_node = " << right_node->export_to_python(out, ctx) << "," << nl;
  }
  ss << ")" << nl << nl;
  if (!str_export) {
    out << ss.str();
    return exported_name;
  }
  else {
    return ss.str();
  }
}

} // namespace API
} // namespace MCell

