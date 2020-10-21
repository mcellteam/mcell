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
#include "gen_wall.h"
#include "../api/wall.h"
#include "../api/geometry_object.h"

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
    throw ValueError("Parameter 'vertices' must be set.");
  }
  if (!is_set(area)) {
    throw ValueError("Parameter 'area' must be set.");
  }
  if (!is_set(normal)) {
    throw ValueError("Parameter 'normal' must be set.");
  }
}

bool GenWall::__eq__(const GenWall& other) const {
  return
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
    wall_index == other.wall_index &&
    vertices == other.vertices &&
    area == other.area &&
    normal == other.normal &&
    is_movable == other.is_movable;
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
  normal = VEC3_UNSET;
  is_movable = true;
}

std::string GenWall::to_str(const std::string ind) const {
  std::stringstream ss;
  ss << get_object_name() << ": " <<
      "\n" << ind + "  " << "geometry_object=" << "(" << ((geometry_object != nullptr) ? geometry_object->to_str(ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "wall_index=" << wall_index << ", " <<
      "vertices=" << vec_nonptr_to_str(vertices, ind + "  ") << ", " <<
      "area=" << area << ", " <<
      "normal=" << normal << ", " <<
      "is_movable=" << is_movable;
  return ss.str();
}

py::class_<Wall> define_pybinding_Wall(py::module& m) {
  return py::class_<Wall, std::shared_ptr<Wall>>(m, "Wall")
      .def(
          py::init<
          >()
      )
      .def("check_semantics", &Wall::check_semantics)
      .def("__str__", &Wall::to_str, py::arg("ind") = std::string(""))
      .def("dump", &Wall::dump)
      .def_property("geometry_object", &Wall::get_geometry_object, &Wall::set_geometry_object)
      .def_property("wall_index", &Wall::get_wall_index, &Wall::set_wall_index)
      .def_property("vertices", &Wall::get_vertices, &Wall::set_vertices)
      .def_property("area", &Wall::get_area, &Wall::set_area)
      .def_property("normal", &Wall::get_normal, &Wall::set_normal)
      .def_property("is_movable", &Wall::get_is_movable, &Wall::set_is_movable)
    ;
}

} // namespace API
} // namespace MCell

