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
  compartment_name = STR_UNSET;
}

std::shared_ptr<ElementaryMolecule> GenElementaryMolecule::copy_elementary_molecule() const {
  if (initialized) {
    throw RuntimeError("Object of class ElementaryMolecule cannot be cloned with 'copy' after this object was used in model initialization.");
  }

  std::shared_ptr<ElementaryMolecule> res = std::make_shared<ElementaryMolecule>(DefaultCtorArgType());
  res->class_name = class_name;
  res->elementary_molecule_type = elementary_molecule_type;
  res->components = components;
  res->compartment_name = compartment_name;

  return res;
}

std::shared_ptr<ElementaryMolecule> GenElementaryMolecule::deepcopy_elementary_molecule(py::dict) const {
  if (initialized) {
    throw RuntimeError("Object of class ElementaryMolecule cannot be cloned with 'deepcopy' after this object was used in model initialization.");
  }

  std::shared_ptr<ElementaryMolecule> res = std::make_shared<ElementaryMolecule>(DefaultCtorArgType());
  res->class_name = class_name;
  res->elementary_molecule_type = is_set(elementary_molecule_type) ? elementary_molecule_type->deepcopy_elementary_molecule_type() : nullptr;
  for (const auto& item: components) {
    res->components.push_back((is_set(item)) ? item->deepcopy_component() : nullptr);
  }
  res->compartment_name = compartment_name;

  return res;
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
    vec_ptr_eq(components, other.components) &&
    compartment_name == other.compartment_name;
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
    true /*components*/ &&
    compartment_name == other.compartment_name;
}

std::string GenElementaryMolecule::to_str(const std::string ind) const {
  std::stringstream ss;
  ss << get_object_name() << ": " <<
      "\n" << ind + "  " << "elementary_molecule_type=" << "(" << ((elementary_molecule_type != nullptr) ? elementary_molecule_type->to_str(ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "components=" << vec_ptr_to_str(components, ind + "  ") << ", " << "\n" << ind + "  " <<
      "compartment_name=" << compartment_name;
  return ss.str();
}

py::class_<ElementaryMolecule> define_pybinding_ElementaryMolecule(py::module& m) {
  return py::class_<ElementaryMolecule, std::shared_ptr<ElementaryMolecule>>(m, "ElementaryMolecule", "Instance of an elementary molecule type. A BNGL complex is composed of elementary molecules.")
      .def(
          py::init<
            std::shared_ptr<ElementaryMoleculeType>,
            const std::vector<std::shared_ptr<Component>>,
            const std::string&
          >(),
          py::arg("elementary_molecule_type"),
          py::arg("components") = std::vector<std::shared_ptr<Component>>(),
          py::arg("compartment_name") = STR_UNSET
      )
      .def("check_semantics", &ElementaryMolecule::check_semantics)
      .def("__copy__", &ElementaryMolecule::copy_elementary_molecule)
      .def("__deepcopy__", &ElementaryMolecule::deepcopy_elementary_molecule, py::arg("memo"))
      .def("__str__", &ElementaryMolecule::to_str, py::arg("ind") = std::string(""))
      .def("__eq__", &ElementaryMolecule::__eq__, py::arg("other"))
      .def("to_bngl_str", &ElementaryMolecule::to_bngl_str, py::arg("with_compartment") = true, "Creates a string that corresponds to its BNGL representation\n- with_compartment: Include compartment name in the returned BNGL string.\n\n")
      .def("dump", &ElementaryMolecule::dump)
      .def_property("elementary_molecule_type", &ElementaryMolecule::get_elementary_molecule_type, &ElementaryMolecule::set_elementary_molecule_type, "Reference to the type of this elementary molecule.")
      .def_property("components", &ElementaryMolecule::get_components, &ElementaryMolecule::set_components, py::return_value_policy::reference, "List of component instances. Not all components need to be specified \nin case when this elementary molecule is used in a pattern.\n")
      .def_property("compartment_name", &ElementaryMolecule::get_compartment_name, &ElementaryMolecule::set_compartment_name, "Optional BNGL compartment name for this elemenrary molecule. If a 2D/surface compartment is specified, the elementary moelcule must be of surface type. If a 3D/volume compartment is specified, the elementary moelcule must be of volume type.")
    ;
}

std::string GenElementaryMolecule::export_to_python(std::ostream& out, PythonExportContext& ctx) {
  if (!export_even_if_already_exported() && ctx.already_exported(this)) {
    return ctx.get_exported_name(this);
  }
  std::string exported_name = "elementary_molecule_" + std::to_string(ctx.postinc_counter("elementary_molecule"));
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
  ss << "m.ElementaryMolecule(" << nl;
  ss << ind << "elementary_molecule_type = " << elementary_molecule_type->export_to_python(out, ctx) << "," << nl;
  if (components != std::vector<std::shared_ptr<Component>>() && !skip_vectors_export()) {
    ss << ind << "components = " << export_vec_components(out, ctx, exported_name) << "," << nl;
  }
  if (compartment_name != STR_UNSET) {
    ss << ind << "compartment_name = " << "'" << compartment_name << "'" << "," << nl;
  }
  ss << ")" << nl << nl;
  if (!str_export) {
    out << ss.str();
    return exported_name;
  }
  else {
    return ss.str();
  }
}

std::string GenElementaryMolecule::export_vec_components(std::ostream& out, PythonExportContext& ctx, const std::string& parent_name) {
  // does not print the array itself to 'out' and returns the whole list
  std::stringstream ss;
  ss << "[";
  for (size_t i = 0; i < components.size(); i++) {
    const auto& item = components[i];
    if (i == 0) {
      ss << " ";
    }
    else if (i % 16 == 0) {
      ss << "\n  ";
    }
    if (!item->skip_python_export()) {
      std::string name = item->export_to_python(out, ctx);
      ss << name << ", ";
    }
  }
  ss << "]";
  return ss.str();
}

} // namespace API
} // namespace MCell

