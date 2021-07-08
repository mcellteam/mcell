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

std::shared_ptr<WallWallHitInfo> GenWallWallHitInfo::copy_wall_wall_hit_info() const {
  std::shared_ptr<WallWallHitInfo> res = std::make_shared<WallWallHitInfo>(DefaultCtorArgType());
  res->class_name = class_name;
  res->wall1 = wall1;
  res->wall2 = wall2;

  return res;
}

std::shared_ptr<WallWallHitInfo> GenWallWallHitInfo::deepcopy_wall_wall_hit_info(py::dict) const {
  std::shared_ptr<WallWallHitInfo> res = std::make_shared<WallWallHitInfo>(DefaultCtorArgType());
  res->class_name = class_name;
  res->wall1 = is_set(wall1) ? wall1->deepcopy_wall() : nullptr;
  res->wall2 = is_set(wall2) ? wall2->deepcopy_wall() : nullptr;

  return res;
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

std::string GenWallWallHitInfo::to_str(const bool all_details, const std::string ind) const {
  std::stringstream ss;
  ss << get_object_name() << ": " <<
      "\n" << ind + "  " << "wall1=" << "(" << ((wall1 != nullptr) ? wall1->to_str(all_details, ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "wall2=" << "(" << ((wall2 != nullptr) ? wall2->to_str(all_details, ind + "  ") : "null" ) << ")";
  return ss.str();
}

py::class_<WallWallHitInfo> define_pybinding_WallWallHitInfo(py::module& m) {
  return py::class_<WallWallHitInfo, std::shared_ptr<WallWallHitInfo>>(m, "WallWallHitInfo", "This class is used in the return type of Model.apply_vertex_moves.\nContains pair of walls that collided.\n")
      .def(
          py::init<
          >()
      )
      .def("check_semantics", &WallWallHitInfo::check_semantics)
      .def("__copy__", &WallWallHitInfo::copy_wall_wall_hit_info)
      .def("__deepcopy__", &WallWallHitInfo::deepcopy_wall_wall_hit_info, py::arg("memo"))
      .def("__str__", &WallWallHitInfo::to_str, py::arg("all_details") = false, py::arg("ind") = std::string(""))
      .def("__eq__", &WallWallHitInfo::__eq__, py::arg("other"))
      .def("dump", &WallWallHitInfo::dump)
      .def_property("wall1", &WallWallHitInfo::get_wall1, &WallWallHitInfo::set_wall1, "First colliding wall.")
      .def_property("wall2", &WallWallHitInfo::get_wall2, &WallWallHitInfo::set_wall2, "Second colliding wall.")
    ;
}

} // namespace API
} // namespace MCell

