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
#include "gen_elementary_molecule.h"
#include "api/elementary_molecule.h"
#include "api/component.h"
#include "api/elementary_molecule_type.h"

namespace MCell {
namespace API {

void GenElementaryMolecule::check_semantics() const {
  if (!is_set(elementary_molecule_type)) {
    throw ValueError("Parameter 'elementary_molecule_type' must be set.");
  }
}

void GenElementaryMolecule::set_initialized() {
  if (is_set(elementary_molecule_type)) {
    elementary_molecule_type->set_initialized();
  }
  vec_set_initialized(components);
  initialized = true;
}

void GenElementaryMolecule::set_all_attributes_as_default_or_unset() {
  class_name = "ElementaryMolecule";
  elementary_molecule_type = nullptr;
  components = std::vector<std::shared_ptr<Component>>();
}

bool GenElementaryMolecule::__eq__(const ElementaryMolecule& other) const {
  return
    (
      (is_set(elementary_molecule_type)) ?
        (is_set(other.elementary_molecule_type) ?
          (elementary_molecule_type->__eq__(*other.elementary_molecule_type)) : 
          false
        ) :
        (is_set(other.elementary_molecule_type) ?
          false :
          true
        )
     )  &&
    vec_ptr_eq(components, other.components);
}

bool GenElementaryMolecule::eq_nonarray_attributes(const ElementaryMolecule& other, const bool ignore_name) const {
  return
    (
      (is_set(elementary_molecule_type)) ?
        (is_set(other.elementary_molecule_type) ?
          (elementary_molecule_type->__eq__(*other.elementary_molecule_type)) : 
          false
        ) :
        (is_set(other.elementary_molecule_type) ?
          false :
          true
        )
     )  &&
    true /*components*/;
}

std::string GenElementaryMolecule::to_str(const std::string ind) const {
  std::stringstream ss;
  ss << get_object_name() << ": " <<
      "\n" << ind + "  " << "elementary_molecule_type=" << "(" << ((elementary_molecule_type != nullptr) ? elementary_molecule_type->to_str(ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "components=" << vec_ptr_to_str(components, ind + "  ");
  return ss.str();
}

py::class_<ElementaryMolecule> define_pybinding_ElementaryMolecule(py::module& m) {
  return py::class_<ElementaryMolecule, std::shared_ptr<ElementaryMolecule>>(m, "ElementaryMolecule")
      .def(
          py::init<
            std::shared_ptr<ElementaryMoleculeType>,
            const std::vector<std::shared_ptr<Component>>
          >(),
          py::arg("elementary_molecule_type"),
          py::arg("components") = std::vector<std::shared_ptr<Component>>()
      )
      .def("check_semantics", &ElementaryMolecule::check_semantics)
      .def("__str__", &ElementaryMolecule::to_str, py::arg("ind") = std::string(""))
      .def("__eq__", &ElementaryMolecule::__eq__, py::arg("other"))
      .def("to_bngl_str", &ElementaryMolecule::to_bngl_str)
      .def("dump", &ElementaryMolecule::dump)
      .def_property("elementary_molecule_type", &ElementaryMolecule::get_elementary_molecule_type, &ElementaryMolecule::set_elementary_molecule_type)
      .def_property("components", &ElementaryMolecule::get_components, &ElementaryMolecule::set_components)
    ;
}

} // namespace API
} // namespace MCell

