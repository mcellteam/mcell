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
#include "api/python_export_utils.h"
#include "gen_geometry_object.h"
#include "api/geometry_object.h"
#include "api/initial_surface_release.h"
#include "api/region.h"
#include "api/surface_class.h"
#include "api/surface_region.h"

namespace MCell {
namespace API {

void GenGeometryObject::check_semantics() const {
  if (!is_set(name)) {
    throw ValueError("Parameter 'name' must be set.");
  }
  if (!is_set(vertex_list)) {
    throw ValueError("Parameter 'vertex_list' must be set and the value must not be an empty list.");
  }
  if (!is_set(wall_list)) {
    throw ValueError("Parameter 'wall_list' must be set and the value must not be an empty list.");
  }
}

void GenGeometryObject::set_initialized() {
  vec_set_initialized(surface_regions);
  if (is_set(surface_class)) {
    surface_class->set_initialized();
  }
  vec_set_initialized(initial_surface_releases);
  if (is_set(left_node)) {
    left_node->set_initialized();
  }
  if (is_set(right_node)) {
    right_node->set_initialized();
  }
  initialized = true;
}

void GenGeometryObject::set_all_attributes_as_default_or_unset() {
  class_name = "GeometryObject";
  name = STR_UNSET;
  vertex_list = std::vector<std::vector<float_t>>();
  wall_list = std::vector<std::vector<int>>();
  is_bngl_compartment = false;
  surface_compartment_name = STR_UNSET;
  surface_regions = std::vector<std::shared_ptr<SurfaceRegion>>();
  surface_class = nullptr;
  initial_surface_releases = std::vector<std::shared_ptr<InitialSurfaceRelease>>();
  node_type = RegionNodeType::UNSET;
  left_node = nullptr;
  right_node = nullptr;
}

bool GenGeometryObject::__eq__(const GeometryObject& other) const {
  return
    name == other.name &&
    vertex_list == other.vertex_list &&
    wall_list == other.wall_list &&
    is_bngl_compartment == other.is_bngl_compartment &&
    surface_compartment_name == other.surface_compartment_name &&
    vec_ptr_eq(surface_regions, other.surface_regions) &&
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

bool GenGeometryObject::eq_nonarray_attributes(const GeometryObject& other, const bool ignore_name) const {
  return
    (ignore_name || name == other.name) &&
    true /*vertex_list*/ &&
    true /*wall_list*/ &&
    is_bngl_compartment == other.is_bngl_compartment &&
    surface_compartment_name == other.surface_compartment_name &&
    true /*surface_regions*/ &&
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

std::string GenGeometryObject::to_str(const std::string ind) const {
  std::stringstream ss;
  ss << get_object_name() << ": " <<
      "name=" << name << ", " <<
      "vertex_list=" << vec_nonptr_to_str(vertex_list, ind + "  ") << ", " <<
      "wall_list=" << vec_nonptr_to_str(wall_list, ind + "  ") << ", " <<
      "is_bngl_compartment=" << is_bngl_compartment << ", " <<
      "surface_compartment_name=" << surface_compartment_name << ", " <<
      "\n" << ind + "  " << "surface_regions=" << vec_ptr_to_str(surface_regions, ind + "  ") << ", " << "\n" << ind + "  " <<
      "surface_class=" << "(" << ((surface_class != nullptr) ? surface_class->to_str(ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "initial_surface_releases=" << vec_ptr_to_str(initial_surface_releases, ind + "  ") << ", " << "\n" << ind + "  " <<
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
            const bool,
            const std::string&,
            const std::vector<std::shared_ptr<SurfaceRegion>>,
            std::shared_ptr<SurfaceClass>,
            const std::vector<std::shared_ptr<InitialSurfaceRelease>>,
            const RegionNodeType,
            std::shared_ptr<Region>,
            std::shared_ptr<Region>
          >(),
          py::arg("name"),
          py::arg("vertex_list"),
          py::arg("wall_list"),
          py::arg("is_bngl_compartment") = false,
          py::arg("surface_compartment_name") = STR_UNSET,
          py::arg("surface_regions") = std::vector<std::shared_ptr<SurfaceRegion>>(),
          py::arg("surface_class") = nullptr,
          py::arg("initial_surface_releases") = std::vector<std::shared_ptr<InitialSurfaceRelease>>(),
          py::arg("node_type") = RegionNodeType::UNSET,
          py::arg("left_node") = nullptr,
          py::arg("right_node") = nullptr
      )
      .def("check_semantics", &GeometryObject::check_semantics)
      .def("__str__", &GeometryObject::to_str, py::arg("ind") = std::string(""))
      .def("__eq__", &GeometryObject::__eq__, py::arg("other"))
      .def("translate", &GeometryObject::translate, py::arg("move"))
      .def("dump", &GeometryObject::dump)
      .def_property("name", &GeometryObject::get_name, &GeometryObject::set_name)
      .def_property("vertex_list", &GeometryObject::get_vertex_list, &GeometryObject::set_vertex_list)
      .def_property("wall_list", &GeometryObject::get_wall_list, &GeometryObject::set_wall_list)
      .def_property("is_bngl_compartment", &GeometryObject::get_is_bngl_compartment, &GeometryObject::set_is_bngl_compartment)
      .def_property("surface_compartment_name", &GeometryObject::get_surface_compartment_name, &GeometryObject::set_surface_compartment_name)
      .def_property("surface_regions", &GeometryObject::get_surface_regions, &GeometryObject::set_surface_regions)
      .def_property("surface_class", &GeometryObject::get_surface_class, &GeometryObject::set_surface_class)
      .def_property("initial_surface_releases", &GeometryObject::get_initial_surface_releases, &GeometryObject::set_initial_surface_releases)
    ;
}

std::string GenGeometryObject::export_to_python(std::ostream& out, PythonExportContext& ctx) {
  if (!export_even_if_already_exported() && ctx.already_exported(this)) {
    return ctx.get_exported_name(this);
  }
  std::string exported_name = std::string("geometry_object") + "_" + (is_set(name) ? fix_id(name) : std::to_string(ctx.postinc_counter("geometry_object")));
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
  ss << "m.GeometryObject(" << nl;
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
  ss << ind << "vertex_list = " << export_vec_vertex_list(out, ctx, exported_name) << "," << nl;
  ss << ind << "wall_list = " << export_vec_wall_list(out, ctx, exported_name) << "," << nl;
  if (is_bngl_compartment != false) {
    ss << ind << "is_bngl_compartment = " << is_bngl_compartment << "," << nl;
  }
  if (surface_compartment_name != STR_UNSET) {
    ss << ind << "surface_compartment_name = " << "'" << surface_compartment_name << "'" << "," << nl;
  }
  if (surface_regions != std::vector<std::shared_ptr<SurfaceRegion>>() && !skip_vectors_export()) {
    ss << ind << "surface_regions = " << export_vec_surface_regions(out, ctx, exported_name) << "," << nl;
  }
  if (is_set(surface_class)) {
    ss << ind << "surface_class = " << surface_class->export_to_python(out, ctx) << "," << nl;
  }
  if (initial_surface_releases != std::vector<std::shared_ptr<InitialSurfaceRelease>>() && !skip_vectors_export()) {
    ss << ind << "initial_surface_releases = " << export_vec_initial_surface_releases(out, ctx, exported_name) << "," << nl;
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

std::string GenGeometryObject::export_vec_vertex_list(std::ostream& out, PythonExportContext& ctx, const std::string& parent_name) {
  // does not print the array itself to 'out' and returns the whole list
  std::stringstream ss;
  ss << "[";
  for (size_t i = 0; i < vertex_list.size(); i++) {
    const auto& item = vertex_list[i];
    if (i == 0) {
      ss << " ";
    }
    else if (i % 16 == 0) {
      ss << "\n  ";
    }
    ss << "[";
    for (const auto& value: item) {
      ss << f_to_str(value) << ", ";
    }
    ss << "], ";
  }
  ss << "]";
  return ss.str();
}

std::string GenGeometryObject::export_vec_wall_list(std::ostream& out, PythonExportContext& ctx, const std::string& parent_name) {
  // does not print the array itself to 'out' and returns the whole list
  std::stringstream ss;
  ss << "[";
  for (size_t i = 0; i < wall_list.size(); i++) {
    const auto& item = wall_list[i];
    if (i == 0) {
      ss << " ";
    }
    else if (i % 16 == 0) {
      ss << "\n  ";
    }
    ss << "[";
    for (const auto& value: item) {
      ss << value << ", ";
    }
    ss << "], ";
  }
  ss << "]";
  return ss.str();
}

std::string GenGeometryObject::export_vec_surface_regions(std::ostream& out, PythonExportContext& ctx, const std::string& parent_name) {
  // does not print the array itself to 'out' and returns the whole list
  std::stringstream ss;
  ss << "[";
  for (size_t i = 0; i < surface_regions.size(); i++) {
    const auto& item = surface_regions[i];
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

std::string GenGeometryObject::export_vec_initial_surface_releases(std::ostream& out, PythonExportContext& ctx, const std::string& parent_name) {
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

