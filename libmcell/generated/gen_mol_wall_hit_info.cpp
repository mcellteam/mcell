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
#include "gen_mol_wall_hit_info.h"
#include "api\mol_wall_hit_info.h"
#include "api\geometry_object.h"

namespace MCell {
namespace API {

bool GenMolWallHitInfo::__eq__(const MolWallHitInfo& other) const {
  return
    molecule_id == other.molecule_id &&
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
    time == other.time &&
    pos3d == other.pos3d &&
    time_before_hit == other.time_before_hit &&
    pos3d_before_hit == other.pos3d_before_hit;
}

bool GenMolWallHitInfo::eq_nonarray_attributes(const MolWallHitInfo& other, const bool ignore_name) const {
  return
    molecule_id == other.molecule_id &&
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
    time == other.time &&
    pos3d == other.pos3d &&
    time_before_hit == other.time_before_hit &&
    pos3d_before_hit == other.pos3d_before_hit;
}

std::string GenMolWallHitInfo::to_str(const std::string ind) const {
  std::stringstream ss;
  ss << "MolWallHitInfo" << ": " <<
      "molecule_id=" << molecule_id << ", " <<
      "\n" << ind + "  " << "geometry_object=" << "(" << ((geometry_object != nullptr) ? geometry_object->to_str(ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "wall_index=" << wall_index << ", " <<
      "time=" << time << ", " <<
      "pos3d=" << pos3d << ", " <<
      "time_before_hit=" << time_before_hit << ", " <<
      "pos3d_before_hit=" << pos3d_before_hit;
  return ss.str();
}

py::class_<MolWallHitInfo> define_pybinding_MolWallHitInfo(py::module& m) {
  return py::class_<MolWallHitInfo, std::shared_ptr<MolWallHitInfo>>(m, "MolWallHitInfo")
      .def(
          py::init<
          >()
      )
      .def("__str__", &MolWallHitInfo::to_str, py::arg("ind") = std::string(""))
      .def("__eq__", &MolWallHitInfo::__eq__, py::arg("other"))
      .def("dump", &MolWallHitInfo::dump)
      .def_property("molecule_id", &MolWallHitInfo::get_molecule_id, &MolWallHitInfo::set_molecule_id)
      .def_property("geometry_object", &MolWallHitInfo::get_geometry_object, &MolWallHitInfo::set_geometry_object)
      .def_property("wall_index", &MolWallHitInfo::get_wall_index, &MolWallHitInfo::set_wall_index)
      .def_property("time", &MolWallHitInfo::get_time, &MolWallHitInfo::set_time)
      .def_property("pos3d", &MolWallHitInfo::get_pos3d, &MolWallHitInfo::set_pos3d)
      .def_property("time_before_hit", &MolWallHitInfo::get_time_before_hit, &MolWallHitInfo::set_time_before_hit)
      .def_property("pos3d_before_hit", &MolWallHitInfo::get_pos3d_before_hit, &MolWallHitInfo::set_pos3d_before_hit)
    ;
}

} // namespace API
} // namespace MCell

