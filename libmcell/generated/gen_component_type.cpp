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
#include "gen_component_type.h"
#include "../api/component_type.h"
#include "../api/component_instance.h"

namespace MCell {
namespace API {

SemRes GenComponentType::check_semantics(std::ostream& out) const{
  if (!is_set(name)) {
    out << get_object_name() << ": Parameter 'name' must be set.\n";
    return SemRes::ERROR;
  }
  return SemRes::OK;
}

std::string GenComponentType::to_str() const{
  std::stringstream ss;
  ss << get_object_name() << ": " <<
      "name=" << name << ", " <<
      "states=" << states;
  return ss.str();
}

py::class_<ComponentType> define_pybinding_ComponentType(py::module& m) {
  return py::class_<ComponentType>(m, "ComponentType")
      .def(
          py::init<
            const std::string&,
            const std::vector<std::string>
          >()
,          py::arg("name"),
          py::arg("states") = std::vector<std::string>()
        )
      .def("check_semantics", &ComponentType::check_semantics_cerr)
      .def("__str__", &ComponentType::to_str)
      .def("inst", &ComponentType::inst, py::arg("state"), py::arg("bond"))
      .def("dump", &ComponentType::dump)
      .def_property("name", &ComponentType::get_name, &ComponentType::set_name)
      .def_property("states", &ComponentType::get_states, &ComponentType::set_states)
    ;
}

} // namespace API
} // namespace MCell

