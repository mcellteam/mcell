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
#include "gen_subsystem.h"
#include "api/subsystem.h"
#include "api/elementary_molecule_type.h"
#include "api/reaction_rule.h"
#include "api/species.h"
#include "api/surface_class.h"

namespace MCell {
namespace API {

bool GenSubsystem::__eq__(const Subsystem& other) const {
  return
    vec_ptr_eq(species, other.species) &&
    vec_ptr_eq(reaction_rules, other.reaction_rules) &&
    vec_ptr_eq(surface_classes, other.surface_classes) &&
    vec_ptr_eq(elementary_molecule_types, other.elementary_molecule_types);
}

bool GenSubsystem::eq_nonarray_attributes(const Subsystem& other, const bool ignore_name) const {
  return
    true /*species*/ &&
    true /*reaction_rules*/ &&
    true /*surface_classes*/ &&
    true /*elementary_molecule_types*/;
}

std::string GenSubsystem::to_str(const std::string ind) const {
  std::stringstream ss;
  ss << "Subsystem" << ": " <<
      "\n" << ind + "  " << "species=" << vec_ptr_to_str(species, ind + "  ") << ", " << "\n" << ind + "  " <<
      "reaction_rules=" << vec_ptr_to_str(reaction_rules, ind + "  ") << ", " << "\n" << ind + "  " <<
      "surface_classes=" << vec_ptr_to_str(surface_classes, ind + "  ") << ", " << "\n" << ind + "  " <<
      "elementary_molecule_types=" << vec_ptr_to_str(elementary_molecule_types, ind + "  ");
  return ss.str();
}

py::class_<Subsystem> define_pybinding_Subsystem(py::module& m) {
  return py::class_<Subsystem, std::shared_ptr<Subsystem>>(m, "Subsystem")
      .def(
          py::init<
            const std::vector<std::shared_ptr<Species>>,
            const std::vector<std::shared_ptr<ReactionRule>>,
            const std::vector<std::shared_ptr<SurfaceClass>>,
            const std::vector<std::shared_ptr<ElementaryMoleculeType>>
          >(),
          py::arg("species") = std::vector<std::shared_ptr<Species>>(),
          py::arg("reaction_rules") = std::vector<std::shared_ptr<ReactionRule>>(),
          py::arg("surface_classes") = std::vector<std::shared_ptr<SurfaceClass>>(),
          py::arg("elementary_molecule_types") = std::vector<std::shared_ptr<ElementaryMoleculeType>>()
      )
      .def("__str__", &Subsystem::to_str, py::arg("ind") = std::string(""))
      .def("__eq__", &Subsystem::__eq__, py::arg("other"))
      .def("add_species", &Subsystem::add_species, py::arg("s"))
      .def("find_species", &Subsystem::find_species, py::arg("name"))
      .def("add_reaction_rule", &Subsystem::add_reaction_rule, py::arg("r"))
      .def("find_reaction_rule", &Subsystem::find_reaction_rule, py::arg("name"))
      .def("add_surface_class", &Subsystem::add_surface_class, py::arg("sc"))
      .def("find_surface_class", &Subsystem::find_surface_class, py::arg("name"))
      .def("add_elementary_molecule_type", &Subsystem::add_elementary_molecule_type, py::arg("mt"))
      .def("find_elementary_molecule_type", &Subsystem::find_elementary_molecule_type, py::arg("name"))
      .def("load_bngl_molecule_types_and_reaction_rules", &Subsystem::load_bngl_molecule_types_and_reaction_rules, py::arg("file_name"), py::arg("parameter_overrides") = std::map<std::string, float_t>())
      .def("dump", &Subsystem::dump)
      .def_property("species", &Subsystem::get_species, &Subsystem::set_species)
      .def_property("reaction_rules", &Subsystem::get_reaction_rules, &Subsystem::set_reaction_rules)
      .def_property("surface_classes", &Subsystem::get_surface_classes, &Subsystem::set_surface_classes)
      .def_property("elementary_molecule_types", &Subsystem::get_elementary_molecule_types, &Subsystem::set_elementary_molecule_types)
    ;
}

std::string GenSubsystem::export_to_python(std::ostream& out, PythonExportContext& ctx) const {
  std::string exported_name = "subsystem";

  std::stringstream ss;
  ss << exported_name << " = m.Subsystem(\n";
  if (species != std::vector<std::shared_ptr<Species>>()) {
    ss << "  species = " << export_vec_species(out, ctx, exported_name) << ",\n";
  }
  if (reaction_rules != std::vector<std::shared_ptr<ReactionRule>>()) {
    ss << "  reaction_rules = " << export_vec_reaction_rules(out, ctx, exported_name) << ",\n";
  }
  if (surface_classes != std::vector<std::shared_ptr<SurfaceClass>>()) {
    ss << "  surface_classes = " << export_vec_surface_classes(out, ctx, exported_name) << ",\n";
  }
  if (elementary_molecule_types != std::vector<std::shared_ptr<ElementaryMoleculeType>>()) {
    ss << "  elementary_molecule_types = " << export_vec_elementary_molecule_types(out, ctx, exported_name) << ",\n";
  }
  ss << ")\n\n";
  out << ss.str();
  return exported_name;
}

std::string GenSubsystem::export_vec_species(std::ostream& out, PythonExportContext& ctx, const std::string& parent_name) const {
  // prints vector into 'out' and returns name of the generated list
  std::stringstream ss;
  std::string exported_name = parent_name + "_species";
  ss << exported_name << " = [\n";
  for (size_t i = 0; i < species.size(); i++) {
    const auto& item = species[i];
    if (i == 0) {
      ss << "  ";
    }
    else if (i % 16 == 0) {
      ss << "\n  ";
    }
    if (!item->skip_python_export()) {
      std::string name = item->export_to_python(out, ctx);
      ss << name << ", ";
    }
  }
  ss << "\n]\n\n";
  out << ss.str();
  return exported_name;
}

std::string GenSubsystem::export_vec_reaction_rules(std::ostream& out, PythonExportContext& ctx, const std::string& parent_name) const {
  // prints vector into 'out' and returns name of the generated list
  std::stringstream ss;
  std::string exported_name = parent_name + "_reaction_rules";
  ss << exported_name << " = [\n";
  for (size_t i = 0; i < reaction_rules.size(); i++) {
    const auto& item = reaction_rules[i];
    if (i == 0) {
      ss << "  ";
    }
    else if (i % 16 == 0) {
      ss << "\n  ";
    }
    if (!item->skip_python_export()) {
      std::string name = item->export_to_python(out, ctx);
      ss << name << ", ";
    }
  }
  ss << "\n]\n\n";
  out << ss.str();
  return exported_name;
}

std::string GenSubsystem::export_vec_surface_classes(std::ostream& out, PythonExportContext& ctx, const std::string& parent_name) const {
  // prints vector into 'out' and returns name of the generated list
  std::stringstream ss;
  std::string exported_name = parent_name + "_surface_classes";
  ss << exported_name << " = [\n";
  for (size_t i = 0; i < surface_classes.size(); i++) {
    const auto& item = surface_classes[i];
    if (i == 0) {
      ss << "  ";
    }
    else if (i % 16 == 0) {
      ss << "\n  ";
    }
    if (!item->skip_python_export()) {
      std::string name = item->export_to_python(out, ctx);
      ss << name << ", ";
    }
  }
  ss << "\n]\n\n";
  out << ss.str();
  return exported_name;
}

std::string GenSubsystem::export_vec_elementary_molecule_types(std::ostream& out, PythonExportContext& ctx, const std::string& parent_name) const {
  // prints vector into 'out' and returns name of the generated list
  std::stringstream ss;
  std::string exported_name = parent_name + "_elementary_molecule_types";
  ss << exported_name << " = [\n";
  for (size_t i = 0; i < elementary_molecule_types.size(); i++) {
    const auto& item = elementary_molecule_types[i];
    if (i == 0) {
      ss << "  ";
    }
    else if (i % 16 == 0) {
      ss << "\n  ";
    }
    if (!item->skip_python_export()) {
      std::string name = item->export_to_python(out, ctx);
      ss << name << ", ";
    }
  }
  ss << "\n]\n\n";
  out << ss.str();
  return exported_name;
}

} // namespace API
} // namespace MCell

