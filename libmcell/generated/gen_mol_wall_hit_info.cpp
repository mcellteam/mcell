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
#include "api/mol_wall_hit_info.h"
#include "api/geometry_object.h"

namespace MCell {
namespace API {

std::shared_ptr<MolWallHitInfo> GenMolWallHitInfo::copy_mol_wall_hit_info() const {
  std::shared_ptr<MolWallHitInfo> res = std::make_shared<MolWallHitInfo>(DefaultCtorArgType());
  res->molecule_id = molecule_id;
  res->geometry_object = geometry_object;
  res->wall_index = wall_index;
  res->time = time;
  res->pos3d = pos3d;
  res->time_before_hit = time_before_hit;
  res->pos3d_before_hit = pos3d_before_hit;

  return res;
}

std::shared_ptr<MolWallHitInfo> GenMolWallHitInfo::deepcopy_mol_wall_hit_info(py::dict) const {
  std::shared_ptr<MolWallHitInfo> res = std::make_shared<MolWallHitInfo>(DefaultCtorArgType());
  res->molecule_id = molecule_id;
  res->geometry_object = is_set(geometry_object) ? geometry_object->deepcopy_geometry_object() : nullptr;
  res->wall_index = wall_index;
  res->time = time;
  res->pos3d = pos3d;
  res->time_before_hit = time_before_hit;
  res->pos3d_before_hit = pos3d_before_hit;

  return res;
}

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

std::string GenMolWallHitInfo::to_str(const bool all_details, const std::string ind) const {
  std::stringstream ss;
  ss << "MolWallHitInfo" << ": " <<
      "molecule_id=" << molecule_id << ", " <<
      "\n" << ind + "  " << "geometry_object=" << "(" << ((geometry_object != nullptr) ? geometry_object->to_str(all_details, ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "wall_index=" << wall_index << ", " <<
      "time=" << time << ", " <<
      "pos3d=" << pos3d << ", " <<
      "time_before_hit=" << time_before_hit << ", " <<
      "pos3d_before_hit=" << pos3d_before_hit;
  return ss.str();
}

py::class_<MolWallHitInfo> define_pybinding_MolWallHitInfo(py::module& m) {
  return py::class_<MolWallHitInfo, std::shared_ptr<MolWallHitInfo>>(m, "MolWallHitInfo", "Data structure passed to a callback function registered through\nModel.register_mol_wall_hit_callback.\n \n")
      .def(
          py::init<
          >()
      )
      .def("__copy__", &MolWallHitInfo::copy_mol_wall_hit_info)
      .def("__deepcopy__", &MolWallHitInfo::deepcopy_mol_wall_hit_info, py::arg("memo"))
      .def("__str__", &MolWallHitInfo::to_str, py::arg("all_details") = false, py::arg("ind") = std::string(""))
      .def("__eq__", &MolWallHitInfo::__eq__, py::arg("other"))
      .def("dump", &MolWallHitInfo::dump)
      .def_property("molecule_id", &MolWallHitInfo::get_molecule_id, &MolWallHitInfo::set_molecule_id, "Id of molecule that hit the wall.")
      .def_property("geometry_object", &MolWallHitInfo::get_geometry_object, &MolWallHitInfo::set_geometry_object, "Object that was hit.")
      .def_property("wall_index", &MolWallHitInfo::get_wall_index, &MolWallHitInfo::set_wall_index, "Index of the wall belonging to the geometry_object.")
      .def_property("time", &MolWallHitInfo::get_time, &MolWallHitInfo::set_time, "Time of the hit.")
      .def_property("pos3d", &MolWallHitInfo::get_pos3d, &MolWallHitInfo::set_pos3d, "Position of the hit.")
      .def_property("time_before_hit", &MolWallHitInfo::get_time_before_hit, &MolWallHitInfo::set_time_before_hit, "The time when the molecule started to diffuse towards the hit wall. \nIt is either the start of the molecule's diffusion or \nwhen the molecule reflected from another wall.\n  \n")
      .def_property("pos3d_before_hit", &MolWallHitInfo::get_pos3d_before_hit, &MolWallHitInfo::set_pos3d_before_hit, "Position of the molecule at time_before_hit.")
    ;
}

} // namespace API
} // namespace MCell

