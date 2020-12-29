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

std::string GenElementaryMolecule::export_to_python(std::ostream& out, PythonExportContext& ctx) const {
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
  ss << ")" << nl << nl;
  if (!str_export) {
    out << ss.str();
    return exported_name;
  }
  else {
    return ss.str();
  }
}

std::string GenElementaryMolecule::export_vec_components(std::ostream& out, PythonExportContext& ctx, const std::string& parent_name) const {
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

