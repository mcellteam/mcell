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
#include "gen_elementary_molecule_instance.h"
#include "api/elementary_molecule_instance.h"
#include "api/component_instance.h"
#include "api/elementary_molecule_type.h"

namespace MCell {
namespace API {

void GenElementaryMoleculeInstance::check_semantics() const {
  if (!is_set(elementary_molecule_type)) {
    throw ValueError("Parameter 'elementary_molecule_type' must be set.");
  }
}

void GenElementaryMoleculeInstance::set_initialized() {
  if (is_set(elementary_molecule_type)) {
    elementary_molecule_type->set_initialized();
  }
  vec_set_initialized(components);
  initialized = true;
}

void GenElementaryMoleculeInstance::set_all_attributes_as_default_or_unset() {
  class_name = "ElementaryMoleculeInstance";
  elementary_molecule_type = nullptr;
  components = std::vector<std::shared_ptr<ComponentInstance>>();
}

bool GenElementaryMoleculeInstance::__eq__(const ElementaryMoleculeInstance& other) const {
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

bool GenElementaryMoleculeInstance::eq_nonarray_attributes(const ElementaryMoleculeInstance& other, const bool ignore_name) const {
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

std::string GenElementaryMoleculeInstance::to_str(const std::string ind) const {
  std::stringstream ss;
  ss << get_object_name() << ": " <<
      "\n" << ind + "  " << "elementary_molecule_type=" << "(" << ((elementary_molecule_type != nullptr) ? elementary_molecule_type->to_str(ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "components=" << vec_ptr_to_str(components, ind + "  ");
  return ss.str();
}

py::class_<ElementaryMoleculeInstance> define_pybinding_ElementaryMoleculeInstance(py::module& m) {
  return py::class_<ElementaryMoleculeInstance, std::shared_ptr<ElementaryMoleculeInstance>>(m, "ElementaryMoleculeInstance")
      .def(
          py::init<
            std::shared_ptr<ElementaryMoleculeType>,
            const std::vector<std::shared_ptr<ComponentInstance>>
          >(),
          py::arg("elementary_molecule_type"),
          py::arg("components") = std::vector<std::shared_ptr<ComponentInstance>>()
      )
      .def("check_semantics", &ElementaryMoleculeInstance::check_semantics)
      .def("__str__", &ElementaryMoleculeInstance::to_str, py::arg("ind") = std::string(""))
      .def("__repr__", &ElementaryMoleculeInstance::to_str, py::arg("ind") = std::string(""))
      .def("__eq__", &ElementaryMoleculeInstance::__eq__, py::arg("other"))
      .def("to_bngl_str", &ElementaryMoleculeInstance::to_bngl_str)
      .def("dump", &ElementaryMoleculeInstance::dump)
      .def_property("elementary_molecule_type", &ElementaryMoleculeInstance::get_elementary_molecule_type, &ElementaryMoleculeInstance::set_elementary_molecule_type)
      .def_property("components", &ElementaryMoleculeInstance::get_components, &ElementaryMoleculeInstance::set_components)
    ;
}

} // namespace API
} // namespace MCell

