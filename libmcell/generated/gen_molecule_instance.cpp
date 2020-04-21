/******************************************************************************
 *
 * Copyright (C) 2020 by
 * The Salk Institute for Biological Studies and
 * Pittsburgh Supercomputing Center, Carnegie Mellon University
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
#include <pybind11/stl.h>
#include "gen_molecule_instance.h"
#include "../api/molecule_instance.h"
#include "../api/component_instance.h"
#include "../api/molecule_type.h"

namespace MCell {
namespace API {

SemRes GenMoleculeInstance::check_semantics(std::ostream& out) const {
  if (!is_set(molecule_type)) {
    out << get_object_name() << ": Parameter 'molecule_type' must be set.\n";
    return SemRes::ERROR;
  }
  return SemRes::OK;
}

std::string GenMoleculeInstance::to_str() const {
  std::stringstream ss;
  ss << get_object_name() << ": " <<
      "molecule_type=" << "(" << ((molecule_type != nullptr) ? molecule_type->to_str() : "null" ) << ")" << ", " <<
      "components=" << vec_ptr_to_str(components);
  return ss.str();
}

py::class_<MoleculeInstance> define_pybinding_MoleculeInstance(py::module& m) {
  return py::class_<MoleculeInstance>(m, "MoleculeInstance")
      .def(
          py::init<
            const MoleculeType*,
            const std::vector<ComponentInstance*>
          >()
,          py::arg("molecule_type"),
          py::arg("components") = std::vector<ComponentInstance*>()
        )
      .def("check_semantics", &MoleculeInstance::check_semantics_cerr)
      .def("__str__", &MoleculeInstance::to_str)
      .def("dump", &MoleculeInstance::dump)
      .def_property("molecule_type", &MoleculeInstance::get_molecule_type, &MoleculeInstance::set_molecule_type)
      .def_property("components", &MoleculeInstance::get_components, &MoleculeInstance::set_components)
    ;
}

} // namespace API
} // namespace MCell

