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
#include "gen_wall_hit_info.h"
#include "../api/wall_hit_info.h"

namespace MCell {
namespace API {

std::string GenWallHitInfo::to_str(const std::string ind) const {
  std::stringstream ss;
  ss << "WallHitInfo" << ": " <<
      "molecule_id=" << molecule_id << ", " <<
      "geometry_object_id=" << geometry_object_id << ", " <<
      "wall_id=" << wall_id << ", " <<
      "time=" << time << ", " <<
      "pos=" << pos << ", " <<
      "time_before_hit=" << time_before_hit << ", " <<
      "pos_before_hit=" << pos_before_hit;
  return ss.str();
}

py::class_<WallHitInfo> define_pybinding_WallHitInfo(py::module& m) {
  return py::class_<WallHitInfo, std::shared_ptr<WallHitInfo>>(m, "WallHitInfo")
      .def(
          py::init<
          >()
      )
      .def("__str__", &WallHitInfo::to_str, py::arg("ind") = std::string(""))
      .def("dump", &WallHitInfo::dump)
      .def_property("molecule_id", &WallHitInfo::get_molecule_id, &WallHitInfo::set_molecule_id)
      .def_property("geometry_object_id", &WallHitInfo::get_geometry_object_id, &WallHitInfo::set_geometry_object_id)
      .def_property("wall_id", &WallHitInfo::get_wall_id, &WallHitInfo::set_wall_id)
      .def_property("time", &WallHitInfo::get_time, &WallHitInfo::set_time)
      .def_property("pos", &WallHitInfo::get_pos, &WallHitInfo::set_pos)
      .def_property("time_before_hit", &WallHitInfo::get_time_before_hit, &WallHitInfo::set_time_before_hit)
      .def_property("pos_before_hit", &WallHitInfo::get_pos_before_hit, &WallHitInfo::set_pos_before_hit)
    ;
}

} // namespace API
} // namespace MCell

