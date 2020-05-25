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
#include "gen_subsystem.h"
#include "../api/subsystem.h"
#include "../api/reaction_rule.h"
#include "../api/species.h"
#include "../api/surface_class.h"

namespace MCell {
namespace API {

std::string GenSubsystem::to_str(const std::string ind) const {
  std::stringstream ss;
  ss << "Subsystem" << ": " <<
      "\n" << ind + "  " << "species=" << vec_ptr_to_str(species, ind + "  ") << ", " << "\n" << ind + "  " <<
      "reaction_rules=" << vec_ptr_to_str(reaction_rules, ind + "  ") << ", " << "\n" << ind + "  " <<
      "surface_classes=" << vec_ptr_to_str(surface_classes, ind + "  ");
  return ss.str();
}

py::class_<Subsystem> define_pybinding_Subsystem(py::module& m) {
  return py::class_<Subsystem, std::shared_ptr<Subsystem>>(m, "Subsystem")
      .def(
          py::init<
          >()
      )
      .def("__str__", &Subsystem::to_str, py::arg("ind") = std::string(""))
      .def("add_species", &Subsystem::add_species, py::arg("s"))
      .def("find_species", &Subsystem::find_species, py::arg("name"))
      .def("add_reaction_rule", &Subsystem::add_reaction_rule, py::arg("r"))
      .def("find_reaction_rule", &Subsystem::find_reaction_rule, py::arg("name"))
      .def("add_surface_class", &Subsystem::add_surface_class, py::arg("sc"))
      .def("dump", &Subsystem::dump)
      .def_property("species", &Subsystem::get_species, &Subsystem::set_species)
      .def_property("reaction_rules", &Subsystem::get_reaction_rules, &Subsystem::set_reaction_rules)
      .def_property("surface_classes", &Subsystem::get_surface_classes, &Subsystem::set_surface_classes)
    ;
}

} // namespace API
} // namespace MCell

