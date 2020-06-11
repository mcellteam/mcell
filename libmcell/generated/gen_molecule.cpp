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
#include "gen_molecule.h"
#include "../api/molecule.h"
#include "../api/species.h"

namespace MCell {
namespace API {

void GenMolecule::check_semantics() const {
}

bool GenMolecule::__eq__(const GenMolecule& other) const {
  return
    name == other.name &&
    id == other.id &&
    (
      (species != nullptr) ?
        ( (other.species != nullptr) ?
          (species->__eq__(*other.species)) : 
          false
        ) :
        ( (other.species != nullptr) ?
          false :
          true
        )
     )  &&
    pos3d == other.pos3d &&
    orientation == other.orientation;
}

void GenMolecule::set_initialized() {
  if (is_set(species)) {
    species->set_initialized();
  }
  initialized = true;
}

std::string GenMolecule::to_str(const std::string ind) const {
  std::stringstream ss;
  ss << get_object_name() << ": " <<
      "id=" << id << ", " <<
      "\n" << ind + "  " << "species=" << "(" << ((species != nullptr) ? species->to_str(ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "pos3d=" << pos3d << ", " <<
      "orientation=" << orientation;
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
      .def("remove", &Molecule::remove)
      .def("dump", &Molecule::dump)
      .def_property("id", &Molecule::get_id, &Molecule::set_id)
      .def_property("species", &Molecule::get_species, &Molecule::set_species)
      .def_property("pos3d", &Molecule::get_pos3d, &Molecule::set_pos3d)
      .def_property("orientation", &Molecule::get_orientation, &Molecule::set_orientation)
    ;
}

} // namespace API
} // namespace MCell

