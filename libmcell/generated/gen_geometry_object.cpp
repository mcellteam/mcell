/******************************************************************************
 *
 * Copyright (C) 2021 by
 * The Salk Institute for Biological Studies
 *
 * Use of this source code is governed by an MIT-style
 * license that can be found in the LICENSE file or at
 * https://opensource.org/licenses/MIT.
 *
******************************************************************************/

#include <sstream>
#include "api/pybind11_stl_include.h"
#include "api/python_export_utils.h"
#include "gen_geometry_object.h"
#include "api/geometry_object.h"
#include "api/color.h"
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

void GenGeometryObject::set_all_attributes_as_default_or_unset() {
  class_name = "GeometryObject";
  name = STR_UNSET;
  vertex_list = std::vector<std::vector<double>>();
  wall_list = std::vector<std::vector<int>>();
  is_bngl_compartment = false;
  surface_compartment_name = STR_UNSET;
  surface_regions = std::vector<std::shared_ptr<SurfaceRegion>>();
  surface_class = nullptr;
  initial_surface_releases = std::vector<std::shared_ptr<InitialSurfaceRelease>>();
  initial_color = nullptr;
  node_type = RegionNodeType::UNSET;
  left_node = nullptr;
  right_node = nullptr;
}

std::shared_ptr<GeometryObject> GenGeometryObject::copy_geometry_object() const {
  std::shared_ptr<GeometryObject> res = std::make_shared<GeometryObject>(DefaultCtorArgType());
  res->class_name = class_name;
  res->name = name;
  res->vertex_list = vertex_list;
  res->wall_list = wall_list;
  res->is_bngl_compartment = is_bngl_compartment;
  res->surface_compartment_name = surface_compartment_name;
  res->surface_regions = surface_regions;
  res->surface_class = surface_class;
  res->initial_surface_releases = initial_surface_releases;
  res->initial_color = initial_color;
  res->node_type = node_type;
  res->left_node = left_node;
  res->right_node = right_node;

  return res;
}

std::shared_ptr<GeometryObject> GenGeometryObject::deepcopy_geometry_object(py::dict) const {
  std::shared_ptr<GeometryObject> res = std::make_shared<GeometryObject>(DefaultCtorArgType());
  res->class_name = class_name;
  res->name = name;
  res->vertex_list = vertex_list;
  res->wall_list = wall_list;
  res->is_bngl_compartment = is_bngl_compartment;
  res->surface_compartment_name = surface_compartment_name;
  for (const auto& item: surface_regions) {
    res->surface_regions.push_back((is_set(item)) ? item->deepcopy_surface_region() : nullptr);
  }
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

std::string GenGeometryObject::to_str(const bool all_details, const std::string ind) const {
  std::stringstream ss;
  ss << get_object_name() << ": " <<
      "name=" << name << ", " <<
      "vertex_list=" << vec_nonptr_to_str(vertex_list, all_details, ind + "  ") << ", " <<
      "wall_list=" << vec_nonptr_to_str(wall_list, all_details, ind + "  ") << ", " <<
      "is_bngl_compartment=" << is_bngl_compartment << ", " <<
      "surface_compartment_name=" << surface_compartment_name << ", " <<
      "\n" << ind + "  " << "surface_regions=" << vec_ptr_to_str(surface_regions, all_details, ind + "  ") << ", " << "\n" << ind + "  " <<
      "surface_class=" << "(" << ((surface_class != nullptr) ? surface_class->to_str(all_details, ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "initial_surface_releases=" << vec_ptr_to_str(initial_surface_releases, all_details, ind + "  ") << ", " << "\n" << ind + "  " <<
      "initial_color=" << "(" << ((initial_color != nullptr) ? initial_color->to_str(all_details, ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "node_type=" << node_type << ", " <<
      "\n" << ind + "  " << "left_node=" << "(" << ((left_node != nullptr) ? left_node->to_str(all_details, ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "right_node=" << "(" << ((right_node != nullptr) ? right_node->to_str(all_details, ind + "  ") : "null" ) << ")";
  return ss.str();
}

py::class_<GeometryObject> define_pybinding_GeometryObject(py::module& m) {
  return py::class_<GeometryObject, Region, std::shared_ptr<GeometryObject>>(m, "GeometryObject", "Class represents geometry objects defined by triangular surface elements.")
      .def(
          py::init<
            const std::string&,
            const std::vector<std::vector<double>>,
            const std::vector<std::vector<int>>,
            const bool,
            const std::string&,
            const std::vector<std::shared_ptr<SurfaceRegion>>,
            std::shared_ptr<SurfaceClass>,
            const std::vector<std::shared_ptr<InitialSurfaceRelease>>,
            std::shared_ptr<Color>,
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
          py::arg("initial_color") = nullptr,
          py::arg("node_type") = RegionNodeType::UNSET,
          py::arg("left_node") = nullptr,
          py::arg("right_node") = nullptr
      )
      .def("check_semantics", &GeometryObject::check_semantics)
      .def("__copy__", &GeometryObject::copy_geometry_object)
      .def("__deepcopy__", &GeometryObject::deepcopy_geometry_object, py::arg("memo"))
      .def("__str__", &GeometryObject::to_str, py::arg("all_details") = false, py::arg("ind") = std::string(""))
      .def("__eq__", &GeometryObject::__eq__, py::arg("other"))
      .def("translate", &GeometryObject::translate, py::arg("move"), "Move object by a specified vector. \nCannot be called after model was initialized.\n\n- move: 3D vector [x, y, z] that will be added to each vertex of this object.\n\n")
      .def("dump", &GeometryObject::dump)
      .def_property("name", &GeometryObject::get_name, &GeometryObject::set_name, "Name of the object. Also represents BNGL compartment name if 'is_bngl_compartment' is True.\n")
      .def_property("vertex_list", &GeometryObject::get_vertex_list, &GeometryObject::set_vertex_list, py::return_value_policy::reference, "List of [x,y,z] triplets specifying positions of individual vertices of each triangle.\n \n")
      .def_property("wall_list", &GeometryObject::get_wall_list, &GeometryObject::set_wall_list, py::return_value_policy::reference, "List of [a,b,c] triplets specifying each wall, individual values are indices into the \nvertex_list attribute.\n")
      .def_property("is_bngl_compartment", &GeometryObject::get_is_bngl_compartment, &GeometryObject::set_is_bngl_compartment, "Set to True if this object represents a 3D BNGL compartment. \nIts name will be then the BNGL compartment name.       \n")
      .def_property("surface_compartment_name", &GeometryObject::get_surface_compartment_name, &GeometryObject::set_surface_compartment_name, "When is_bngl_compartment is True, this attribute can be set to specify its \nmembrane (2D) compartment name.\n")
      .def_property("surface_regions", &GeometryObject::get_surface_regions, &GeometryObject::set_surface_regions, py::return_value_policy::reference, "All surface regions associated with this geometry object.\n")
      .def_property("surface_class", &GeometryObject::get_surface_class, &GeometryObject::set_surface_class, "Surface class for the whole object's surface. It is applied to the whole surface of this object \nexcept for those surface regions that have their specific surface class set explicitly.\n")
      .def_property("initial_surface_releases", &GeometryObject::get_initial_surface_releases, &GeometryObject::set_initial_surface_releases, py::return_value_policy::reference, "Each item in this list defines either density or number of molecules to be released on this surface \nregions when simulation starts.\n")
      .def_property("initial_color", &GeometryObject::get_initial_color, &GeometryObject::set_initial_color, "Initial color for this geometry object. If a surface region has its color set, its value \nis used for the walls of that surface region.\n")
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

