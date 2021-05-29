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
#include "gen_component_type.h"
#include "api/component_type.h"
#include "api/component.h"

namespace MCell {
namespace API {

void GenComponentType::check_semantics() const {
  if (!is_set(name)) {
    throw ValueError("Parameter 'name' must be set.");
  }
}

void GenComponentType::set_initialized() {
  initialized = true;
}

void GenComponentType::set_all_attributes_as_default_or_unset() {
  class_name = "ComponentType";
  name = STR_UNSET;
  states = std::vector<std::string>();
}

std::shared_ptr<ComponentType> GenComponentType::copy_component_type() const {
  if (initialized) {
    throw RuntimeError("Object of class ComponentType cannot be cloned with 'copy' after this object was used in model initialization.");
  }

  std::shared_ptr<ComponentType> res = std::make_shared<ComponentType>(DefaultCtorArgType());
  res->class_name = class_name;
  res->name = name;
  res->states = states;

  return res;
}

std::shared_ptr<ComponentType> GenComponentType::deepcopy_component_type(py::dict) const {
  if (initialized) {
    throw RuntimeError("Object of class ComponentType cannot be cloned with 'deepcopy' after this object was used in model initialization.");
  }

  std::shared_ptr<ComponentType> res = std::make_shared<ComponentType>(DefaultCtorArgType());
  res->class_name = class_name;
  res->name = name;
  res->states = states;

  return res;
}

bool GenComponentType::__eq__(const ComponentType& other) const {
  return
    name == other.name &&
    states == other.states;
}

bool GenComponentType::eq_nonarray_attributes(const ComponentType& other, const bool ignore_name) const {
  return
    (ignore_name || name == other.name) &&
    true /*states*/;
}

std::string GenComponentType::to_str(const std::string ind) const {
  std::stringstream ss;
  ss << get_object_name() << ": " <<
      "name=" << name << ", " <<
      "states=" << vec_nonptr_to_str(states, ind + "  ");
  return ss.str();
}

py::class_<ComponentType> define_pybinding_ComponentType(py::module& m) {
  return py::class_<ComponentType, std::shared_ptr<ComponentType>>(m, "ComponentType", "Multiple functional attributes for each molecule type are described using components. And this class defines a type of a component. For example, proteins have multiple functional substructures such as domains, motifs, and binding sites. These components can be unchanging (called stateless) or exist in one of many different internal states For example, certain binding motifs may have different behaviors depending on whether they are unphosphorylated or phosphorylated.")
      .def(
          py::init<
            const std::string&,
            const std::vector<std::string>
          >(),
          py::arg("name"),
          py::arg("states") = std::vector<std::string>()
      )
      .def("check_semantics", &ComponentType::check_semantics)
      .def("__copy__", &ComponentType::copy_component_type)
      .def("__deepcopy__", &ComponentType::deepcopy_component_type, py::arg("memo"))
      .def("__str__", &ComponentType::to_str, py::arg("ind") = std::string(""))
      .def("__eq__", &ComponentType::__eq__, py::arg("other"))
      .def("inst", py::overload_cast<const std::string&, const int>(&ComponentType::inst), py::arg("state") = "STATE_UNSET", py::arg("bond") = BOND_UNBOUND, "Instantiate a component from this component type.\n- state: Selected state, must be from the list of the allowed states.\n\n- bond: Bond information for the created component instance.\n\n")
      .def("inst", py::overload_cast<const int, const int>(&ComponentType::inst), py::arg("state") = STATE_UNSET_INT, py::arg("bond") = BOND_UNBOUND, "Instantiate a component from this component type.\n- state: Selected state, must be from the list of the allowed, converted to string.\n\n- bond: Bond information for the created component instance.\n\n")
      .def("to_bngl_str", &ComponentType::to_bngl_str, "Creates a string that corresponds to its BNGL representation.")
      .def("dump", &ComponentType::dump)
      .def_property("name", &ComponentType::get_name, &ComponentType::set_name, "Name of this component type.")
      .def_property("states", &ComponentType::get_states, &ComponentType::set_states, py::return_value_policy::reference, "List of states allowed by this component.")
    ;
}

std::string GenComponentType::export_to_python(std::ostream& out, PythonExportContext& ctx) {
  if (!export_even_if_already_exported() && ctx.already_exported(this)) {
    return ctx.get_exported_name(this);
  }
  std::string exported_name = std::string("component_type") + "_" + (is_set(name) ? fix_id(name) : std::to_string(ctx.postinc_counter("component_type")));
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
  ss << "m.ComponentType(" << nl;
  ss << ind << "name = " << "'" << name << "'" << "," << nl;
  if (states != std::vector<std::string>() && !skip_vectors_export()) {
    ss << ind << "states = " << export_vec_states(out, ctx, exported_name) << "," << nl;
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

std::string GenComponentType::export_vec_states(std::ostream& out, PythonExportContext& ctx, const std::string& parent_name) {
  // does not print the array itself to 'out' and returns the whole list
  std::stringstream ss;
  ss << "[";
  for (size_t i = 0; i < states.size(); i++) {
    const auto& item = states[i];
    if (i == 0) {
      ss << " ";
    }
    else if (i % 16 == 0) {
      ss << "\n  ";
    }
    ss << "'" << item << "', ";
  }
  ss << "]";
  return ss.str();
}

} // namespace API
} // namespace MCell

