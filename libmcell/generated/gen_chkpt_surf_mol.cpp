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
#include "gen_chkpt_surf_mol.h"
#include "api/chkpt_surf_mol.h"
#include "api/geometry_object.h"
#include "api/species.h"

namespace MCell {
namespace API {

void GenChkptSurfMol::check_semantics() const {
  if (!is_set(pos)) {
    throw ValueError("Parameter 'pos' must be set.");
  }
  if (!is_set(orientation)) {
    throw ValueError("Parameter 'orientation' must be set.");
  }
  if (!is_set(geometry_object)) {
    throw ValueError("Parameter 'geometry_object' must be set.");
  }
  if (!is_set(wall_index)) {
    throw ValueError("Parameter 'wall_index' must be set.");
  }
  if (!is_set(grid_tile_index)) {
    throw ValueError("Parameter 'grid_tile_index' must be set.");
  }
  if (!is_set(id)) {
    throw ValueError("Parameter 'id' must be set.");
  }
  if (!is_set(species)) {
    throw ValueError("Parameter 'species' must be set.");
  }
  if (!is_set(diffusion_time)) {
    throw ValueError("Parameter 'diffusion_time' must be set.");
  }
  if (!is_set(unimol_rx_time)) {
    throw ValueError("Parameter 'unimol_rx_time' must be set.");
  }
  if (!is_set(birthday)) {
    throw ValueError("Parameter 'birthday' must be set.");
  }
}

void GenChkptSurfMol::set_initialized() {
  if (is_set(geometry_object)) {
    geometry_object->set_initialized();
  }
  if (is_set(species)) {
    species->set_initialized();
  }
  initialized = true;
}

void GenChkptSurfMol::set_all_attributes_as_default_or_unset() {
  class_name = "ChkptSurfMol";
  pos = VEC2_UNSET;
  orientation = Orientation::NOT_SET;
  geometry_object = nullptr;
  wall_index = INT_UNSET;
  grid_tile_index = INT_UNSET;
  id = INT_UNSET;
  species = nullptr;
  diffusion_time = FLT_UNSET;
  unimol_rx_time = FLT_UNSET;
  birthday = FLT_UNSET;
}

bool GenChkptSurfMol::__eq__(const ChkptSurfMol& other) const {
  return
    pos == other.pos &&
    orientation == other.orientation &&
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
    grid_tile_index == other.grid_tile_index &&
    id == other.id &&
    (
      (is_set(species)) ?
        (is_set(other.species) ?
          (species->__eq__(*other.species)) : 
          false
        ) :
        (is_set(other.species) ?
          false :
          true
        )
     )  &&
    diffusion_time == other.diffusion_time &&
    unimol_rx_time == other.unimol_rx_time &&
    birthday == other.birthday;
}

bool GenChkptSurfMol::eq_nonarray_attributes(const ChkptSurfMol& other, const bool ignore_name) const {
  return
    pos == other.pos &&
    orientation == other.orientation &&
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
    grid_tile_index == other.grid_tile_index &&
    id == other.id &&
    (
      (is_set(species)) ?
        (is_set(other.species) ?
          (species->__eq__(*other.species)) : 
          false
        ) :
        (is_set(other.species) ?
          false :
          true
        )
     )  &&
    diffusion_time == other.diffusion_time &&
    unimol_rx_time == other.unimol_rx_time &&
    birthday == other.birthday;
}

std::string GenChkptSurfMol::to_str(const std::string ind) const {
  std::stringstream ss;
  ss << get_object_name() << ": " <<
      "pos=" << pos << ", " <<
      "orientation=" << orientation << ", " <<
      "\n" << ind + "  " << "geometry_object=" << "(" << ((geometry_object != nullptr) ? geometry_object->to_str(ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "wall_index=" << wall_index << ", " <<
      "grid_tile_index=" << grid_tile_index << ", " <<
      "id=" << id << ", " <<
      "\n" << ind + "  " << "species=" << "(" << ((species != nullptr) ? species->to_str(ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "diffusion_time=" << diffusion_time << ", " <<
      "unimol_rx_time=" << unimol_rx_time << ", " <<
      "birthday=" << birthday;
  return ss.str();
}

py::class_<ChkptSurfMol> define_pybinding_ChkptSurfMol(py::module& m) {
  return py::class_<ChkptSurfMol, BaseChkptMol, std::shared_ptr<ChkptSurfMol>>(m, "ChkptSurfMol")
      .def(
          py::init<
            const Vec2&,
            const Orientation,
            std::shared_ptr<GeometryObject>,
            const int,
            const int,
            const int,
            std::shared_ptr<Species>,
            const float_t,
            const float_t,
            const float_t
          >(),
          py::arg("pos"),
          py::arg("orientation"),
          py::arg("geometry_object"),
          py::arg("wall_index"),
          py::arg("grid_tile_index"),
          py::arg("id"),
          py::arg("species"),
          py::arg("diffusion_time"),
          py::arg("unimol_rx_time"),
          py::arg("birthday")
      )
      .def("check_semantics", &ChkptSurfMol::check_semantics)
      .def("__str__", &ChkptSurfMol::to_str, py::arg("ind") = std::string(""))
      .def("__eq__", &ChkptSurfMol::__eq__, py::arg("other"))
      .def("dump", &ChkptSurfMol::dump)
      .def_property("pos", &ChkptSurfMol::get_pos, &ChkptSurfMol::set_pos)
      .def_property("orientation", &ChkptSurfMol::get_orientation, &ChkptSurfMol::set_orientation)
      .def_property("geometry_object", &ChkptSurfMol::get_geometry_object, &ChkptSurfMol::set_geometry_object)
      .def_property("wall_index", &ChkptSurfMol::get_wall_index, &ChkptSurfMol::set_wall_index)
      .def_property("grid_tile_index", &ChkptSurfMol::get_grid_tile_index, &ChkptSurfMol::set_grid_tile_index)
    ;
}

std::string GenChkptSurfMol::export_to_python(std::ostream& out, PythonExportContext& ctx) const {
  if (!export_even_if_already_exported() && ctx.already_exported(this)) {
    return ctx.get_exported_name(this);
  }
  std::string exported_name = "chkpt_surf_mol_" + std::to_string(ctx.postinc_counter("chkpt_surf_mol"));
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
  ss << "m.ChkptSurfMol(" << nl;
  ss << ind << "pos = " << pos << "," << nl;
  ss << ind << "orientation = " << orientation << "," << nl;
  ss << ind << "geometry_object = " << geometry_object->export_to_python(out, ctx) << "," << nl;
  ss << ind << "wall_index = " << wall_index << "," << nl;
  ss << ind << "grid_tile_index = " << grid_tile_index << "," << nl;
  ss << ind << "id = " << id << "," << nl;
  ss << ind << "species = " << species->export_to_python(out, ctx) << "," << nl;
  ss << ind << "diffusion_time = " << diffusion_time << "," << nl;
  ss << ind << "unimol_rx_time = " << unimol_rx_time << "," << nl;
  ss << ind << "birthday = " << birthday << "," << nl;
  ss << ")" << nl << nl;
  if (!str_export) {
    out << ss.str();
    return exported_name;
  }
  else {
    return ss.str();
  }
}

} // namespace API
} // namespace MCell

