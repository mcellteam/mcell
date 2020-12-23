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
  return py::class_<ComponentType, std::shared_ptr<ComponentType>>(m, "ComponentType")
      .def(
          py::init<
            const std::string&,
            const std::vector<std::string>
          >(),
          py::arg("name"),
          py::arg("states") = std::vector<std::string>()
      )
      .def("check_semantics", &ComponentType::check_semantics)
      .def("__str__", &ComponentType::to_str, py::arg("ind") = std::string(""))
      .def("__eq__", &ComponentType::__eq__, py::arg("other"))
      .def("inst", py::overload_cast<const std::string&, const int>(&ComponentType::inst), py::arg("state") = "STATE_UNSET", py::arg("bond") = BOND_UNBOUND)
      .def("inst", py::overload_cast<const int, const int>(&ComponentType::inst), py::arg("state") = STATE_UNSET_INT, py::arg("bond") = BOND_UNBOUND)
      .def("to_bngl_str", &ComponentType::to_bngl_str)
      .def("dump", &ComponentType::dump)
      .def_property("name", &ComponentType::get_name, &ComponentType::set_name)
      .def_property("states", &ComponentType::get_states, &ComponentType::set_states)
    ;
}

std::string GenComponentType::export_to_python(std::ostream& out, PythonExportContext& ctx) const {
  if (ctx.already_exported(this)) {
    return ctx.get_exported_name(this);
  }
  std::string exported_name = fix_id(name);
  ctx.add_exported(this, exported_name);

  std::stringstream ss;
  ss << exported_name << " = ComponentType(\n";
  ss << "  name = " << name << ",\n";
  if (states != std::vector<std::string>()) {
    ss << "  states = " << export_vec_states(out, ctx, exported_name) << ",\n";
  }
  ss << ")\n\n";
  out << ss.str();
  return exported_name;
}

std::string GenComponentType::export_vec_states(std::ostream& out, PythonExportContext& ctx, const std::string& parent_name) const {
  std::string exported_name = parent_name + "_states";
  std::stringstream ss;
  ss << exported_name << " = [\n";
  for (size_t i = 0; i < states.size(); i++) {
    const auto& item = states[i];
    if (i == 0) {
      ss << "  ";
    }
    else if (i % 16 == 0) {
      ss << "\n  ";
    }
    ss << item << ", ";
  }
  ss << "]\n\n";
  out << ss.str();
  return exported_name;
}

} // namespace API
} // namespace MCell

