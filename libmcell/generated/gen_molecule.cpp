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
#include "gen_molecule.h"
#include "api/molecule.h"
#include "api/geometry_object.h"

namespace MCell {
namespace API {

void GenMolecule::check_semantics() const {
}

void GenMolecule::set_initialized() {
  if (is_set(geometry_object)) {
    geometry_object->set_initialized();
  }
  initialized = true;
}

void GenMolecule::set_all_attributes_as_default_or_unset() {
  class_name = "Molecule";
  id = ID_INVALID;
  type = MoleculeType::UNSET;
  species_id = ID_INVALID;
  pos3d = VEC3_UNSET;
  orientation = Orientation::NOT_SET;
  pos2d = VEC2_UNSET;
  geometry_object = nullptr;
  wall_index = -1;
}

bool GenMolecule::__eq__(const Molecule& other) const {
  return
    id == other.id &&
    type == other.type &&
    species_id == other.species_id &&
    pos3d == other.pos3d &&
    orientation == other.orientation &&
    pos2d == other.pos2d &&
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
    wall_index == other.wall_index;
}

bool GenMolecule::eq_nonarray_attributes(const Molecule& other, const bool ignore_name) const {
  return
    id == other.id &&
    type == other.type &&
    species_id == other.species_id &&
    pos3d == other.pos3d &&
    orientation == other.orientation &&
    pos2d == other.pos2d &&
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
    wall_index == other.wall_index;
}

std::string GenMolecule::to_str(const std::string ind) const {
  std::stringstream ss;
  ss << get_object_name() << ": " <<
      "id=" << id << ", " <<
      "type=" << type << ", " <<
      "species_id=" << species_id << ", " <<
      "pos3d=" << pos3d << ", " <<
      "orientation=" << orientation << ", " <<
      "pos2d=" << pos2d << ", " <<
      "\n" << ind + "  " << "geometry_object=" << "(" << ((geometry_object != nullptr) ? geometry_object->to_str(ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "wall_index=" << wall_index;
  return ss.str();
}

py::class_<Molecule> define_pybinding_Molecule(py::module& m) {
  return py::class_<Molecule, std::shared_ptr<Molecule>>(m, "Molecule")
      .def(
          py::init<
          >()
      )
      .def("check_semantics", &Molecule::check_semantics)
      .def("__str__", &Molecule::to_str, py::arg("ind") = std::string(""))
      .def("__eq__", &Molecule::__eq__, py::arg("other"))
      .def("remove", &Molecule::remove)
      .def("dump", &Molecule::dump)
      .def_property("id", &Molecule::get_id, &Molecule::set_id)
      .def_property("type", &Molecule::get_type, &Molecule::set_type)
      .def_property("species_id", &Molecule::get_species_id, &Molecule::set_species_id)
      .def_property("pos3d", &Molecule::get_pos3d, &Molecule::set_pos3d)
      .def_property("orientation", &Molecule::get_orientation, &Molecule::set_orientation)
      .def_property("pos2d", &Molecule::get_pos2d, &Molecule::set_pos2d)
      .def_property("geometry_object", &Molecule::get_geometry_object, &Molecule::set_geometry_object)
      .def_property("wall_index", &Molecule::get_wall_index, &Molecule::set_wall_index)
    ;
}

} // namespace API
} // namespace MCell

