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
#include <pybind11/stl.h>
#include "gen_component_instance.h"
#include "../api/component_instance.h"
#include "../api/component_type.h"

namespace MCell {
namespace API {

void GenComponentInstance::check_semantics() const {
  if (!is_set(component_type)) {
    throw ValueError("Parameter 'component_type' must be set.");
  }
}

std::string GenComponentInstance::to_str(const std::string ind) const {
  std::stringstream ss;
  ss << get_object_name() << ": " <<
      "\n" << ind + "  " << "component_type=" << "(" << ((component_type != nullptr) ? component_type->to_str(ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "state=" << state << ", " <<
      "bond=" << bond;
  return ss.str();
}

py::class_<ComponentInstance> define_pybinding_ComponentInstance(py::module& m) {
  return py::class_<ComponentInstance, std::shared_ptr<ComponentInstance>>(m, "ComponentInstance")
      .def(
          py::init<
            std::shared_ptr<ComponentType>,
            const std::string&,
            const int
          >(),
          py::arg("component_type"),
          py::arg("state") = STATE_UNSET,
          py::arg("bond") = BOND_UNBOUND
      )
      .def("check_semantics", &ComponentInstance::check_semantics)
      .def("__str__", &ComponentInstance::to_str, py::arg("ind") = std::string(""))
      .def("dump", &ComponentInstance::dump)
      .def_property("component_type", &ComponentInstance::get_component_type, &ComponentInstance::set_component_type)
      .def_property("state", &ComponentInstance::get_state, &ComponentInstance::set_state)
      .def_property("bond", &ComponentInstance::get_bond, &ComponentInstance::set_bond)
    ;
}

} // namespace API
} // namespace MCell

