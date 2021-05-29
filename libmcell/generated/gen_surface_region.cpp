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
#include "gen_surface_region.h"
#include "api/surface_region.h"
#include "api/color.h"
#include "api/initial_surface_release.h"
#include "api/region.h"
#include "api/surface_class.h"

namespace MCell {
namespace API {

void GenSurfaceRegion::check_semantics() const {
  if (!is_set(name)) {
    throw ValueError("Parameter 'name' must be set.");
  }
  if (!is_set(wall_indices)) {
    throw ValueError("Parameter 'wall_indices' must be set and the value must not be an empty list.");
  }
}

void GenSurfaceRegion::set_initialized() {
  if (is_set(surface_class)) {
    surface_class->set_initialized();
  }
  vec_set_initialized(initial_surface_releases);
  if (is_set(initial_color)) {
    initial_color->set_initialized();
  }
  if (is_set(left_node)) {
    left_node->set_initialized();
  }
  if (is_set(right_node)) {
    right_node->set_initialized();
  }
  initialized = true;
}

void GenSurfaceRegion::set_all_attributes_as_default_or_unset() {
  class_name = "SurfaceRegion";
  name = STR_UNSET;
  wall_indices = std::vector<int>();
  surface_class = nullptr;
  initial_surface_releases = std::vector<std::shared_ptr<InitialSurfaceRelease>>();
  initial_color = nullptr;
  node_type = RegionNodeType::UNSET;
  left_node = nullptr;
  right_node = nullptr;
}

std::shared_ptr<SurfaceRegion> GenSurfaceRegion::copy_surface_region() const {
  std::shared_ptr<SurfaceRegion> res = std::make_shared<SurfaceRegion>(DefaultCtorArgType());
  res->class_name = class_name;
  res->name = name;
  res->wall_indices = wall_indices;
  res->surface_class = surface_class;
  res->initial_surface_releases = initial_surface_releases;
  res->initial_color = initial_color;
  res->node_type = node_type;
  res->left_node = left_node;
  res->right_node = right_node;

  return res;
}

std::shared_ptr<SurfaceRegion> GenSurfaceRegion::deepcopy_surface_region(py::dict) const {
  std::shared_ptr<SurfaceRegion> res = std::make_shared<SurfaceRegion>(DefaultCtorArgType());
  res->class_name = class_name;
  res->name = name;
  res->wall_indices = wall_indices;
  res->surface_class = is_set(surface_class) ? surface_class->deepcopy_surface_class() : nullptr;
  for (const auto& item: initial_surface_releases) {
    res->initial_surface_releases.push_back((is_set(item)) ? item->deepcopy_initial_surface_release() : nullptr);
  }
  res->initial_color = is_set(initial_color) ? initial_color->deepcopy_color() : nullptr;
  res->node_type = node_type;
  res->left_node = is_set(left_node) ? left_node->deepcopy_region() : nullptr;
  res->right_node = is_set(right_node) ? right_node->deepcopy_region() : nullptr;

  return res;
}

bool GenSurfaceRegion::__eq__(const SurfaceRegion& other) const {
  return
    name == other.name &&
    wall_indices == other.wall_indices &&
    (
      (is_set(surface_class)) ?
        (is_set(other.surface_class) ?
          (surface_class->__eq__(*other.surface_class)) : 
          false
        ) :
        (is_set(other.surface_class) ?
          false :
          true
        )
     )  &&
    vec_ptr_eq(initial_surface_releases, other.initial_surface_releases) &&
    (
      (is_set(initial_color)) ?
        (is_set(other.initial_color) ?
          (initial_color->__eq__(*other.initial_color)) : 
          false
        ) :
        (is_set(other.initial_color) ?
          false :
          true
        )
     )  &&
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

bool GenSurfaceRegion::eq_nonarray_attributes(const SurfaceRegion& other, const bool ignore_name) const {
  return
    (ignore_name || name == other.name) &&
    true /*wall_indices*/ &&
    (
      (is_set(surface_class)) ?
        (is_set(other.surface_class) ?
          (surface_class->__eq__(*other.surface_class)) : 
          false
        ) :
        (is_set(other.surface_class) ?
          false :
          true
        )
     )  &&
    true /*initial_surface_releases*/ &&
    (
      (is_set(initial_color)) ?
        (is_set(other.initial_color) ?
          (initial_color->__eq__(*other.initial_color)) : 
          false
        ) :
        (is_set(other.initial_color) ?
          false :
          true
        )
     )  &&
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

std::string GenSurfaceRegion::to_str(const bool all_details, const std::string ind) const {
  std::stringstream ss;
  ss << get_object_name() << ": " <<
      "name=" << name << ", " <<
      "wall_indices=" << vec_nonptr_to_str(wall_indices, all_details, ind + "  ") << ", " <<
      "\n" << ind + "  " << "surface_class=" << "(" << ((surface_class != nullptr) ? surface_class->to_str(all_details, ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "initial_surface_releases=" << vec_ptr_to_str(initial_surface_releases, all_details, ind + "  ") << ", " << "\n" << ind + "  " <<
      "initial_color=" << "(" << ((initial_color != nullptr) ? initial_color->to_str(all_details, ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "node_type=" << node_type << ", " <<
      "\n" << ind + "  " << "left_node=" << "(" << ((left_node != nullptr) ? left_node->to_str(all_details, ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "right_node=" << "(" << ((right_node != nullptr) ? right_node->to_str(all_details, ind + "  ") : "null" ) << ")";
  return ss.str();
}

py::class_<SurfaceRegion> define_pybinding_SurfaceRegion(py::module& m) {
  return py::class_<SurfaceRegion, Region, std::shared_ptr<SurfaceRegion>>(m, "SurfaceRegion", "Defines a region on the object. The extent of a region is given by the wall_indices list. \nMolecules can be added and surface properties can be set with the optional regional surface commands. \nYou can have an arbitrary number of regions on an object, and they may overlap if\nyou wish. Molecules added to overlapping regions accumulate. Triangles belonging to \nmultiple regions inherit all parent regions’ surface properties. Users\nhave to make sure that in case of overlapped regions their surface properties\nare compatible. \n")
      .def(
          py::init<
            const std::string&,
            const std::vector<int>,
            std::shared_ptr<SurfaceClass>,
            const std::vector<std::shared_ptr<InitialSurfaceRelease>>,
            std::shared_ptr<Color>,
            const RegionNodeType,
            std::shared_ptr<Region>,
            std::shared_ptr<Region>
          >(),
          py::arg("name"),
          py::arg("wall_indices"),
          py::arg("surface_class") = nullptr,
          py::arg("initial_surface_releases") = std::vector<std::shared_ptr<InitialSurfaceRelease>>(),
          py::arg("initial_color") = nullptr,
          py::arg("node_type") = RegionNodeType::UNSET,
          py::arg("left_node") = nullptr,
          py::arg("right_node") = nullptr
      )
      .def("check_semantics", &SurfaceRegion::check_semantics)
      .def("__copy__", &SurfaceRegion::copy_surface_region)
      .def("__deepcopy__", &SurfaceRegion::deepcopy_surface_region, py::arg("memo"))
      .def("__str__", &SurfaceRegion::to_str, py::arg("all_details") = false, py::arg("ind") = std::string(""))
      .def("__eq__", &SurfaceRegion::__eq__, py::arg("other"))
      .def("dump", &SurfaceRegion::dump)
      .def_property("name", &SurfaceRegion::get_name, &SurfaceRegion::set_name, "Name of this region.")
      .def_property("wall_indices", &SurfaceRegion::get_wall_indices, &SurfaceRegion::set_wall_indices, py::return_value_policy::reference, "Surface region must be a part of a GeometryObject, items in this list are indices to \nits wall_list array.\n")
      .def_property("surface_class", &SurfaceRegion::get_surface_class, &SurfaceRegion::set_surface_class, "Optional surface class assigned to this surface region.\nIf not set, it is inherited from the parent geometry object's surface_class.\n")
      .def_property("initial_surface_releases", &SurfaceRegion::get_initial_surface_releases, &SurfaceRegion::set_initial_surface_releases, py::return_value_policy::reference, "Each item of this list defines either density or number of molecules to be released on this surface \nregions when simulation starts.\n")
      .def_property("initial_color", &SurfaceRegion::get_initial_color, &SurfaceRegion::set_initial_color, "Initial color for this specific surface region. If not set, color of the parent's GeometryObject is used.")
    ;
}

std::string GenSurfaceRegion::export_to_python(std::ostream& out, PythonExportContext& ctx) {
  if (!export_even_if_already_exported() && ctx.already_exported(this)) {
    return ctx.get_exported_name(this);
  }
  std::string exported_name = std::string("surface_region") + "_" + (is_set(name) ? fix_id(name) : std::to_string(ctx.postinc_counter("surface_region")));
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
  ss << "m.SurfaceRegion(" << nl;
  if (node_type != RegionNodeType::UNSET) {
    ss << ind << "node_type = " << node_type << "," << nl;
  }
  if (is_set(left_node)) {
    ss << ind << "left_node = " << left_node->export_to_python(out, ctx) << "," << nl;
  }
  if (is_set(right_node)) {
    ss << ind << "right_node = " << right_node->export_to_python(out, ctx) << "," << nl;
  }
  ss << ind << "name = " << "'" << name << "'" << "," << nl;
  ss << ind << "wall_indices = " << export_vec_wall_indices(out, ctx, exported_name) << "," << nl;
  if (is_set(surface_class)) {
    ss << ind << "surface_class = " << surface_class->export_to_python(out, ctx) << "," << nl;
  }
  if (initial_surface_releases != std::vector<std::shared_ptr<InitialSurfaceRelease>>() && !skip_vectors_export()) {
    ss << ind << "initial_surface_releases = " << export_vec_initial_surface_releases(out, ctx, exported_name) << "," << nl;
  }
  if (is_set(initial_color)) {
    ss << ind << "initial_color = " << initial_color->export_to_python(out, ctx) << "," << nl;
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

std::string GenSurfaceRegion::export_vec_wall_indices(std::ostream& out, PythonExportContext& ctx, const std::string& parent_name) {
  // does not print the array itself to 'out' and returns the whole list
  std::stringstream ss;
  ss << "[";
  for (size_t i = 0; i < wall_indices.size(); i++) {
    const auto& item = wall_indices[i];
    if (i == 0) {
      ss << " ";
    }
    else if (i % 16 == 0) {
      ss << "\n  ";
    }
    ss << item << ", ";
  }
  ss << "]";
  return ss.str();
}

std::string GenSurfaceRegion::export_vec_initial_surface_releases(std::ostream& out, PythonExportContext& ctx, const std::string& parent_name) {
  // does not print the array itself to 'out' and returns the whole list
  std::stringstream ss;
  ss << "[";
  for (size_t i = 0; i < initial_surface_releases.size(); i++) {
    const auto& item = initial_surface_releases[i];
    if (i == 0) {
      ss << " ";
    }
    else if (i % 16 == 0) {
      ss << "\n  ";
    }
    if (!item->skip_python_export()) {
      std::string name = item->export_to_python(out, ctx);
      ss << name << ", ";
    }
  }
  ss << "]";
  return ss.str();
}

} // namespace API
} // namespace MCell

