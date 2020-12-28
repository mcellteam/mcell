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
#include "gen_component.h"
#include "api/component.h"
#include "api/component_type.h"

namespace MCell {
namespace API {

void GenComponent::check_semantics() const {
  if (!is_set(component_type)) {
    throw ValueError("Parameter 'component_type' must be set.");
  }
}

void GenComponent::set_initialized() {
  if (is_set(component_type)) {
    component_type->set_initialized();
  }
  initialized = true;
}

void GenComponent::set_all_attributes_as_default_or_unset() {
  class_name = "Component";
  component_type = nullptr;
  state = "STATE_UNSET";
  bond = BOND_UNBOUND;
}

bool GenComponent::__eq__(const Component& other) const {
  return
    (
      (is_set(component_type)) ?
        (is_set(other.component_type) ?
          (component_type->__eq__(*other.component_type)) : 
          false
        ) :
        (is_set(other.component_type) ?
          false :
          true
        )
     )  &&
    state == other.state &&
    bond == other.bond;
}

bool GenComponent::eq_nonarray_attributes(const Component& other, const bool ignore_name) const {
  return
    (
      (is_set(component_type)) ?
        (is_set(other.component_type) ?
          (component_type->__eq__(*other.component_type)) : 
          false
        ) :
        (is_set(other.component_type) ?
          false :
          true
        )
     )  &&
    state == other.state &&
    bond == other.bond;
}

std::string GenComponent::to_str(const std::string ind) const {
  std::stringstream ss;
  ss << get_object_name() << ": " <<
      "\n" << ind + "  " << "component_type=" << "(" << ((component_type != nullptr) ? component_type->to_str(ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "state=" << state << ", " <<
      "bond=" << bond;
  return ss.str();
}

py::class_<Component> define_pybinding_Component(py::module& m) {
  return py::class_<Component, std::shared_ptr<Component>>(m, "Component")
      .def(
          py::init<
            std::shared_ptr<ComponentType>,
            const std::string&,
            const int
          >(),
          py::arg("component_type"),
          py::arg("state") = "STATE_UNSET",
          py::arg("bond") = BOND_UNBOUND
      )
      .def("check_semantics", &Component::check_semantics)
      .def("__str__", &Component::to_str, py::arg("ind") = std::string(""))
      .def("__eq__", &Component::__eq__, py::arg("other"))
      .def("to_bngl_str", &Component::to_bngl_str)
      .def("dump", &Component::dump)
      .def_property("component_type", &Component::get_component_type, &Component::set_component_type)
      .def_property("state", &Component::get_state, &Component::set_state)
      .def_property("bond", &Component::get_bond, &Component::set_bond)
    ;
}

std::string GenComponent::export_to_python(std::ostream& out, PythonExportContext& ctx) const {
  if (ctx.already_exported(this)) {
    return ctx.get_exported_name(this);
  }
  std::string exported_name = "component_" + std::to_string(ctx.postinc_counter("component"));
  ctx.add_exported(this, exported_name);

  bool str_export = export_as_string_without_newlines();
  std::string nl = "";
  std::string ind = " ";
  std::stringstream ss;
  if (!str_export) {
    nl = "\n";
    ind = "    ";
    ss << exported_name << " = ";
  }
  ss << "m.Component(" << nl;
  ss << ind << "component_type = " << component_type->export_to_python(out, ctx) << "," << nl;
  if (state != "STATE_UNSET") {
    ss << ind << "state = " << "'" << name << "'" << "," << nl;
  }
  if (bond != BOND_UNBOUND) {
    ss << ind << "bond = " << bond << "," << nl;
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

} // namespace API
} // namespace MCell

