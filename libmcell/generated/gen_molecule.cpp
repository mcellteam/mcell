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

std::shared_ptr<Molecule> GenMolecule::copy_molecule() const {
  std::shared_ptr<Molecule> res = std::make_shared<Molecule>(DefaultCtorArgType());
  res->class_name = class_name;
  res->id = id;
  res->type = type;
  res->species_id = species_id;
  res->pos3d = pos3d;
  res->orientation = orientation;
  res->pos2d = pos2d;
  res->geometry_object = geometry_object;
  res->wall_index = wall_index;

  return res;
}

std::shared_ptr<Molecule> GenMolecule::deepcopy_molecule(py::dict) const {
  std::shared_ptr<Molecule> res = std::make_shared<Molecule>(DefaultCtorArgType());
  res->class_name = class_name;
  res->id = id;
  res->type = type;
  res->species_id = species_id;
  res->pos3d = pos3d;
  res->orientation = orientation;
  res->pos2d = pos2d;
  res->geometry_object = is_set(geometry_object) ? geometry_object->deepcopy_geometry_object() : nullptr;
  res->wall_index = wall_index;

  return res;
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

std::string GenMolecule::to_str(const bool all_details, const std::string ind) const {
  std::stringstream ss;
  ss << get_object_name() << ": " <<
      "id=" << id << ", " <<
      "type=" << type << ", " <<
      "species_id=" << species_id << ", " <<
      "pos3d=" << pos3d << ", " <<
      "orientation=" << orientation << ", " <<
      "pos2d=" << pos2d << ", " <<
      "\n" << ind + "  " << "geometry_object=" << "(" << ((geometry_object != nullptr) ? geometry_object->to_str(all_details, ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "wall_index=" << wall_index;
  return ss.str();
}

py::class_<Molecule> define_pybinding_Molecule(py::module& m) {
  return py::class_<Molecule, std::shared_ptr<Molecule>>(m, "Molecule", "Representation of a molecule obtained from Model \nduring simulation obtained through Model.get_molecule.\nChanges through changing attributes of this object are not allowed except \nfor complete removal of this molecule.   \n")
      .def(
          py::init<
          >()
      )
      .def("check_semantics", &Molecule::check_semantics)
      .def("__copy__", &Molecule::copy_molecule)
      .def("__deepcopy__", &Molecule::deepcopy_molecule, py::arg("memo"))
      .def("__str__", &Molecule::to_str, py::arg("all_details") = false, py::arg("ind") = std::string(""))
      .def("__eq__", &Molecule::__eq__, py::arg("other"))
      .def("remove", &Molecule::remove, "Removes this molecule from simulation. Any subsequent modifications\nof this molecules won't have any effect.\n")
      .def("dump", &Molecule::dump)
      .def_property("id", &Molecule::get_id, &Molecule::set_id, "Unique id of this molecule. MCell assigns this unique id to each created \nmolecule. All reactions change ID of molecules even in reactions such as \nA@CP -> A@EC.\n")
      .def_property("type", &Molecule::get_type, &Molecule::set_type, "Type of this molecule, either volume or surface. \n")
      .def_property("species_id", &Molecule::get_species_id, &Molecule::set_species_id, "Species id of this molecule.\nThe species id value is only temporary and can be invalidated by simulating an \niteration.\n")
      .def_property("pos3d", &Molecule::get_pos3d, &Molecule::set_pos3d, "Contains position of a molecule in 3D space.        \n")
      .def_property("orientation", &Molecule::get_orientation, &Molecule::set_orientation, "Contains orientation for surface molecule. Volume molecules \nhave always orientation set to Orientation.NONE.\n")
      .def_property("pos2d", &Molecule::get_pos2d, &Molecule::set_pos2d, "Set only for surface molecules. Position on a wall in UV coordinates \nrelative to the triangle of the wall.\n        \n")
      .def_property("geometry_object", &Molecule::get_geometry_object, &Molecule::set_geometry_object, "Set only for surface molecules.\nObject on whose surface is the molecule located.\n")
      .def_property("wall_index", &Molecule::get_wall_index, &Molecule::set_wall_index, "Set only for surface molecules.\nIndex of wall belonging to the geometry_object where is the \nmolecule located. \n   \n")
    ;
}

} // namespace API
} // namespace MCell

