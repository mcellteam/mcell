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
#include "gen_wall_wall_hit_info.h"
#include "api/wall_wall_hit_info.h"
#include "api/wall.h"

namespace MCell {
namespace API {

void GenWallWallHitInfo::check_semantics() const {
  if (!is_set(wall1)) {
    throw ValueError("Parameter 'wall1' must be set.");
  }
  if (!is_set(wall2)) {
    throw ValueError("Parameter 'wall2' must be set.");
  }
}

void GenWallWallHitInfo::set_initialized() {
  if (is_set(wall1)) {
    wall1->set_initialized();
  }
  if (is_set(wall2)) {
    wall2->set_initialized();
  }
  initialized = true;
}

void GenWallWallHitInfo::set_all_attributes_as_default_or_unset() {
  class_name = "WallWallHitInfo";
  wall1 = nullptr;
  wall2 = nullptr;
}

bool GenWallWallHitInfo::__eq__(const WallWallHitInfo& other) const {
  return
    (
      (is_set(wall1)) ?
        (is_set(other.wall1) ?
          (wall1->__eq__(*other.wall1)) : 
          false
        ) :
        (is_set(other.wall1) ?
          false :
          true
        )
     )  &&
    (
      (is_set(wall2)) ?
        (is_set(other.wall2) ?
          (wall2->__eq__(*other.wall2)) : 
          false
        ) :
        (is_set(other.wall2) ?
          false :
          true
        )
     ) ;
}

bool GenWallWallHitInfo::eq_nonarray_attributes(const WallWallHitInfo& other, const bool ignore_name) const {
  return
    (
      (is_set(wall1)) ?
        (is_set(other.wall1) ?
          (wall1->__eq__(*other.wall1)) : 
          false
        ) :
        (is_set(other.wall1) ?
          false :
          true
        )
     )  &&
    (
      (is_set(wall2)) ?
        (is_set(other.wall2) ?
          (wall2->__eq__(*other.wall2)) : 
          false
        ) :
        (is_set(other.wall2) ?
          false :
          true
        )
     ) ;
}

std::string GenWallWallHitInfo::to_str(const std::string ind) const {
  std::stringstream ss;
  ss << get_object_name() << ": " <<
      "\n" << ind + "  " << "wall1=" << "(" << ((wall1 != nullptr) ? wall1->to_str(ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "wall2=" << "(" << ((wall2 != nullptr) ? wall2->to_str(ind + "  ") : "null" ) << ")";
  return ss.str();
}

py::class_<WallWallHitInfo> define_pybinding_WallWallHitInfo(py::module& m) {
  return py::class_<WallWallHitInfo, std::shared_ptr<WallWallHitInfo>>(m, "WallWallHitInfo", "This class is used in the return type of Model.apply_vertex_moves.\nContains pair of walls that collided.\n")
      .def(
          py::init<
          >()
      )
      .def("check_semantics", &WallWallHitInfo::check_semantics)
      .def("__str__", &WallWallHitInfo::to_str, py::arg("ind") = std::string(""))
      .def("__eq__", &WallWallHitInfo::__eq__, py::arg("other"))
      .def("dump", &WallWallHitInfo::dump)
      .def_property("wall1", &WallWallHitInfo::get_wall1, &WallWallHitInfo::set_wall1, "First colliding wall.")
      .def_property("wall2", &WallWallHitInfo::get_wall2, &WallWallHitInfo::set_wall2, "Second colliding wall.")
    ;
}

} // namespace API
} // namespace MCell

