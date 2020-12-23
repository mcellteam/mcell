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
#include "api/python_export.h"
#include "gen_surface_region.h"
#include "api/surface_region.h"
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
  node_type = RegionNodeType::UNSET;
  left_node = nullptr;
  right_node = nullptr;
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

std::string GenSurfaceRegion::to_str(const std::string ind) const {
  std::stringstream ss;
  ss << get_object_name() << ": " <<
      "name=" << name << ", " <<
      "wall_indices=" << vec_nonptr_to_str(wall_indices, ind + "  ") << ", " <<
      "\n" << ind + "  " << "surface_class=" << "(" << ((surface_class != nullptr) ? surface_class->to_str(ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "initial_surface_releases=" << vec_ptr_to_str(initial_surface_releases, ind + "  ") << ", " << "\n" << ind + "  " <<
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
            std::shared_ptr<SurfaceClass>,
            const std::vector<std::shared_ptr<InitialSurfaceRelease>>,
            const RegionNodeType,
            std::shared_ptr<Region>,
            std::shared_ptr<Region>
          >(),
          py::arg("name"),
          py::arg("wall_indices"),
          py::arg("surface_class") = nullptr,
          py::arg("initial_surface_releases") = std::vector<std::shared_ptr<InitialSurfaceRelease>>(),
          py::arg("node_type") = RegionNodeType::UNSET,
          py::arg("left_node") = nullptr,
          py::arg("right_node") = nullptr
      )
      .def("check_semantics", &SurfaceRegion::check_semantics)
      .def("__str__", &SurfaceRegion::to_str, py::arg("ind") = std::string(""))
      .def("__eq__", &SurfaceRegion::__eq__, py::arg("other"))
      .def("dump", &SurfaceRegion::dump)
      .def_property("name", &SurfaceRegion::get_name, &SurfaceRegion::set_name)
      .def_property("wall_indices", &SurfaceRegion::get_wall_indices, &SurfaceRegion::set_wall_indices)
      .def_property("surface_class", &SurfaceRegion::get_surface_class, &SurfaceRegion::set_surface_class)
      .def_property("initial_surface_releases", &SurfaceRegion::get_initial_surface_releases, &SurfaceRegion::set_initial_surface_releases)
    ;
}

std::string GenSurfaceRegion::export_to_python(std::ostream& out, PythonExportContext& ctx) const {
  if (ctx.already_exported(this)) {
    return ctx.get_exported_name(this);
  }
  std::string exported_name = fix_id(name);
  ctx.add_exported(this, exported_name);

  std::stringstream ss;
  ss << exported_name << " = SurfaceRegion(\n";
  ss << "  name = " << name << ",\n";
  ss << "  wall_indices = " << export_vec_wall_indices(out, ctx) << ",\n";
  if (is_set(surface_class)) {
    ss << "  surface_class = " << surface_class->export_to_python(out, ctx) << ",\n";
  }
  if (initial_surface_releases != std::vector<std::shared_ptr<InitialSurfaceRelease>>()) {
    ss << "  initial_surface_releases = " << export_vec_initial_surface_releases(out, ctx) << ",\n";
  }
  if (node_type != RegionNodeType::UNSET) {
    ss << "  node_type = " << node_type << ",\n";
  }
  if (is_set(left_node)) {
    ss << "  left_node = " << left_node->export_to_python(out, ctx) << ",\n";
  }
  if (is_set(right_node)) {
    ss << "  right_node = " << right_node->export_to_python(out, ctx) << ",\n";
  }
  ss << ")\n\n";
  out << ss.str();
  return exported_name;
}

std::string GenSurfaceRegion::export_vec_wall_indices(std::ostream& out, PythonExportContext& ctx) const {
  return ""; //TODO
}

std::string GenSurfaceRegion::export_vec_initial_surface_releases(std::ostream& out, PythonExportContext& ctx) const {
  return ""; //TODO
}

} // namespace API
} // namespace MCell

