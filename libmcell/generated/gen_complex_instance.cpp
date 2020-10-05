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
#include "gen_complex_instance.h"
#include "../api/complex_instance.h"
#include "../api/elementary_molecule_instance.h"
#include "../api/species.h"

namespace MCell {
namespace API {

void GenComplexInstance::check_semantics() const {
}

bool GenComplexInstance::__eq__(const GenComplexInstance& other) const {
  return
    name == other.name &&
    bngl_string == other.bngl_string &&
    vec_ptr_eq(elementary_molecule_instances, other.elementary_molecule_instances) &&
    orientation == other.orientation;
}

void GenComplexInstance::set_initialized() {
  vec_set_initialized(elementary_molecule_instances);
  initialized = true;
}

void GenComplexInstance::set_all_attributes_as_default_or_unset() {
  class_name = "ComplexInstance";
  bngl_string = STR_UNSET;
  elementary_molecule_instances = std::vector<std::shared_ptr<ElementaryMoleculeInstance>>();
  orientation = Orientation::NONE;
}

std::string GenComplexInstance::to_str(const std::string ind) const {
  std::stringstream ss;
  ss << get_object_name() << ": " <<
      "bngl_string=" << bngl_string << ", " <<
      "\n" << ind + "  " << "elementary_molecule_instances=" << vec_ptr_to_str(elementary_molecule_instances, ind + "  ") << ", " << "\n" << ind + "  " <<
      "orientation=" << orientation;
  return ss.str();
}

py::class_<ComplexInstance> define_pybinding_ComplexInstance(py::module& m) {
  return py::class_<ComplexInstance, std::shared_ptr<ComplexInstance>>(m, "ComplexInstance")
      .def(
          py::init<
            const std::string&,
            const std::vector<std::shared_ptr<ElementaryMoleculeInstance>>,
            const Orientation
          >(),
          py::arg("bngl_string") = STR_UNSET,
          py::arg("elementary_molecule_instances") = std::vector<std::shared_ptr<ElementaryMoleculeInstance>>(),
          py::arg("orientation") = Orientation::NONE
      )
      .def("check_semantics", &ComplexInstance::check_semantics)
      .def("__str__", &ComplexInstance::to_str, py::arg("ind") = std::string(""))
      .def("to_bngl_str", &ComplexInstance::to_bngl_str)
      .def("as_species", &ComplexInstance::as_species)
      .def("dump", &ComplexInstance::dump)
      .def_property("bngl_string", &ComplexInstance::get_bngl_string, &ComplexInstance::set_bngl_string)
      .def_property("elementary_molecule_instances", &ComplexInstance::get_elementary_molecule_instances, &ComplexInstance::set_elementary_molecule_instances)
      .def_property("orientation", &ComplexInstance::get_orientation, &ComplexInstance::set_orientation)
    ;
}

} // namespace API
} // namespace MCell

