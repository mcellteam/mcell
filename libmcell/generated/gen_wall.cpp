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
#include "gen_wall.h"
#include "api/wall.h"
#include "api/geometry_object.h"

namespace MCell {
namespace API {

void GenWall::check_semantics() const {
  if (!is_set(geometry_object)) {
    throw ValueError("Parameter 'geometry_object' must be set.");
  }
  if (!is_set(wall_index)) {
    throw ValueError("Parameter 'wall_index' must be set.");
  }
  if (!is_set(vertices)) {
    throw ValueError("Parameter 'vertices' must be set and the value must not be an empty list.");
  }
  if (!is_set(area)) {
    throw ValueError("Parameter 'area' must be set.");
  }
  if (!is_set(unit_normal)) {
    throw ValueError("Parameter 'unit_normal' must be set.");
  }
}

void GenWall::set_initialized() {
  if (is_set(geometry_object)) {
    geometry_object->set_initialized();
  }
  initialized = true;
}

void GenWall::set_all_attributes_as_default_or_unset() {
  class_name = "Wall";
  geometry_object = nullptr;
  wall_index = INT_UNSET;
  vertices = std::vector<Vec3>();
  area = FLT_UNSET;
  unit_normal = VEC3_UNSET;
  is_movable = true;
}

std::shared_ptr<Wall> GenWall::copy_wall() const {
  std::shared_ptr<Wall> res = std::make_shared<Wall>(DefaultCtorArgType());
  res->class_name = class_name;
  res->geometry_object = geometry_object;
  res->wall_index = wall_index;
  res->vertices = vertices;
  res->area = area;
  res->unit_normal = unit_normal;
  res->is_movable = is_movable;

  return res;
}

std::shared_ptr<Wall> GenWall::deepcopy_wall(py::dict) const {
  std::shared_ptr<Wall> res = std::make_shared<Wall>(DefaultCtorArgType());
  res->class_name = class_name;
  res->geometry_object = is_set(geometry_object) ? geometry_object->deepcopy_geometry_object() : nullptr;
  res->wall_index = wall_index;
  res->vertices = vertices;
  res->area = area;
  res->unit_normal = unit_normal;
  res->is_movable = is_movable;

  return res;
}

bool GenWall::__eq__(const Wall& other) const {
  return
    (
      (is_set(geometry_object)) ?
        (is_set(other.geometry_object) ?
          (geometry_object->__eq__(*other.geometry_object)) : 
          false
        ) :
        (is_set(other.geometry_object) ?
          false :
          true
        )
     )  &&
    wall_index == other.wall_index &&
    vertices == other.vertices &&
    area == other.area &&
    unit_normal == other.unit_normal &&
    is_movable == other.is_movable;
}

bool GenWall::eq_nonarray_attributes(const Wall& other, const bool ignore_name) const {
  return
    (
      (is_set(geometry_object)) ?
        (is_set(other.geometry_object) ?
          (geometry_object->__eq__(*other.geometry_object)) : 
          false
        ) :
        (is_set(other.geometry_object) ?
          false :
          true
        )
     )  &&
    wall_index == other.wall_index &&
    true /*vertices*/ &&
    area == other.area &&
    unit_normal == other.unit_normal &&
    is_movable == other.is_movable;
}

std::string GenWall::to_str(const bool all_details, const std::string ind) const {
  std::stringstream ss;
  ss << get_object_name() << ": " <<
      "\n" << ind + "  " << "geometry_object=" << "(" << ((geometry_object != nullptr) ? geometry_object->to_str(all_details, ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "wall_index=" << wall_index << ", " <<
      "vertices=" << vec_nonptr_to_str(vertices, all_details, ind + "  ") << ", " <<
      "area=" << area << ", " <<
      "unit_normal=" << unit_normal << ", " <<
      "is_movable=" << is_movable;
  return ss.str();
}

py::class_<Wall> define_pybinding_Wall(py::module& m) {
  return py::class_<Wall, std::shared_ptr<Wall>>(m, "Wall", "Constant representation of wall of a geometry object.\nChanges through changing attributes of this object are not allowed\nexcept for the attribute is_movable.\n")
      .def(
          py::init<
          >()
      )
      .def("check_semantics", &Wall::check_semantics)
      .def("__copy__", &Wall::copy_wall)
      .def("__deepcopy__", &Wall::deepcopy_wall, py::arg("memo"))
      .def("__str__", &Wall::to_str, py::arg("all_details") = false, py::arg("ind") = std::string(""))
      .def("__eq__", &Wall::__eq__, py::arg("other"))
      .def("dump", &Wall::dump)
      .def_property("geometry_object", &Wall::get_geometry_object, &Wall::set_geometry_object, "Object to which this wall belongs.")
      .def_property("wall_index", &Wall::get_wall_index, &Wall::set_wall_index, "Index of this wall in the object to which this wall belongs.")
      .def_property("vertices", &Wall::get_vertices, &Wall::set_vertices, py::return_value_policy::reference, "Vertices of the triangle that represents this wall.")
      .def_property("area", &Wall::get_area, &Wall::set_area, "Area of the wall in um^2.")
      .def_property("unit_normal", &Wall::get_unit_normal, &Wall::set_unit_normal, "Normal of this wall with unit length of 1 um.\nThere is also a method Model.get_wall_unit_normal that allows to \nretrieve just the normal value without the need to prepare this \nwhole Wall object.  \n")
      .def_property("is_movable", &Wall::get_is_movable, &Wall::set_is_movable, "If True, whis wall can be moved through Model.apply_vertex_moves,\nif False, wall moves are ignored. \nCan be set during simulation.\n")
    ;
}

} // namespace API
} // namespace MCell

