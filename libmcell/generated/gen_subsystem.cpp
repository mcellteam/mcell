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
  return py::class_<Subsystem, std::shared_ptr<Subsystem>>(m, "Subsystem", "Subsystem usually defines a reaction network. It is a collection of \nspecies and reaction rules that use these species. \nThe main motivation for introducing such an object to MCell4 is to have \na class independent on that particular initial model state and observables that \nonly contains reactions. This way, one can define independent reusable subsystems\nand possibly merge them together when creating a model that includes multiple reaction \nnetworks. \n")
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
      .def("add_species", &Subsystem::add_species, py::arg("s"), "Add a reference to a Species object to the species list.\n- s\n")
      .def("find_species", &Subsystem::find_species, py::arg("name"), "Find a Species object using name in the species list. \nReturns None if no such species is found.\n\n- name\n")
      .def("add_reaction_rule", &Subsystem::add_reaction_rule, py::arg("r"), "Add a reference to a ReactionRule object to the reaction_rules list.\n- r\n")
      .def("find_reaction_rule", &Subsystem::find_reaction_rule, py::arg("name"), "Find a ReactionRule object using name in the reaction_rules list. \nReturns None if no such reaction rule is found.\n\n- name\n")
      .def("add_surface_class", &Subsystem::add_surface_class, py::arg("sc"), "Add a reference to a SurfaceClass object to the surface_classes list.\n- sc\n")
      .def("find_surface_class", &Subsystem::find_surface_class, py::arg("name"), "Find a SurfaceClass object using name in the surface_classes list. \nReturns None if no such surface class is found.\n\n- name\n")
      .def("add_elementary_molecule_type", &Subsystem::add_elementary_molecule_type, py::arg("mt"), "Add a reference to an ElementaryMoleculeType object to the elementary_molecule_types list.\n- mt\n")
      .def("find_elementary_molecule_type", &Subsystem::find_elementary_molecule_type, py::arg("name"), "Find an ElementaryMoleculeType object using name in the elementary_molecule_types list. \nReturns None if no such elementary molecule type is found.\n\n- name\n")
      .def("load_bngl_molecule_types_and_reaction_rules", &Subsystem::load_bngl_molecule_types_and_reaction_rules, py::arg("file_name"), py::arg("parameter_overrides") = std::map<std::string, float_t>(), "Parses a BNGL file, only reads molecule types and reaction rules sections, \ni.e. ignores observables and seed species. \nParameter values are evaluated and the result value is directly used.  \nCompartments names are stored in rxn rules as strings because compartments belong \nto geometry objects and the subsystem is independent on specific geometry.\nHowever, the compartments and their objects must be defined before initialization.\n\n- file_name: Path to the BNGL file to be loaded.\n\n- parameter_overrides: For each key k in the parameter_overrides, if it is defined in the BNGL's parameters section,\nits value is ignored and instead value parameter_overrides[k] is used.\n\n\n")
      .def("dump", &Subsystem::dump)
      .def_property("species", &Subsystem::get_species, &Subsystem::set_species, "List of species to be included in the model for initialization.\nUsed usually only for simple species (species that are defined using a\nsingle molecule type without components such as 'A').\nOther species may be created inside simulation  \n")
      .def_property("reaction_rules", &Subsystem::get_reaction_rules, &Subsystem::set_reaction_rules)
      .def_property("surface_classes", &Subsystem::get_surface_classes, &Subsystem::set_surface_classes)
      .def_property("elementary_molecule_types", &Subsystem::get_elementary_molecule_types, &Subsystem::set_elementary_molecule_types, "Contains list of elementary molecule types with their diffusion constants and other information. \nPopulated when a BNGL file is loaded and also on initialization from Species objects present in \nthe species list.\n")
    ;
}

std::string GenSubsystem::export_to_python(std::ostream& out, PythonExportContext& ctx) {
  std::string exported_name = "subsystem";

  bool str_export = export_as_string_without_newlines();
  std::string nl = "";
  std::string ind = " ";
  std::stringstream ss;
  if (!str_export) {
    nl = "\n";
    ind = "    ";
    ss << exported_name << " = ";
  }
  ss << "m.Subsystem(" << nl;
  if (species != std::vector<std::shared_ptr<Species>>() && !skip_vectors_export()) {
    ss << ind << "species = " << export_vec_species(out, ctx, exported_name) << "," << nl;
  }
  if (reaction_rules != std::vector<std::shared_ptr<ReactionRule>>() && !skip_vectors_export()) {
    ss << ind << "reaction_rules = " << export_vec_reaction_rules(out, ctx, exported_name) << "," << nl;
  }
  if (surface_classes != std::vector<std::shared_ptr<SurfaceClass>>() && !skip_vectors_export()) {
    ss << ind << "surface_classes = " << export_vec_surface_classes(out, ctx, exported_name) << "," << nl;
  }
  if (elementary_molecule_types != std::vector<std::shared_ptr<ElementaryMoleculeType>>() && !skip_vectors_export()) {
    ss << ind << "elementary_molecule_types = " << export_vec_elementary_molecule_types(out, ctx, exported_name) << "," << nl;
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

std::string GenSubsystem::export_vec_species(std::ostream& out, PythonExportContext& ctx, const std::string& parent_name) {
  // prints vector into 'out' and returns name of the generated list
  std::stringstream ss;
  std::string exported_name;
  if (parent_name != ""){
    exported_name = parent_name+ "_species";
  }
  else {
    exported_name = "species";
  }

  ss << exported_name << " = [\n";
  for (size_t i = 0; i < species.size(); i++) {
    const auto& item = species[i];
    if (i == 0) {
      ss << "    ";
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

std::string GenSubsystem::export_vec_reaction_rules(std::ostream& out, PythonExportContext& ctx, const std::string& parent_name) {
  // prints vector into 'out' and returns name of the generated list
  std::stringstream ss;
  std::string exported_name;
  if (parent_name != ""){
    exported_name = parent_name+ "_reaction_rules";
  }
  else {
    exported_name = "reaction_rules";
  }

  ss << exported_name << " = [\n";
  for (size_t i = 0; i < reaction_rules.size(); i++) {
    const auto& item = reaction_rules[i];
    if (i == 0) {
      ss << "    ";
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

std::string GenSubsystem::export_vec_surface_classes(std::ostream& out, PythonExportContext& ctx, const std::string& parent_name) {
  // prints vector into 'out' and returns name of the generated list
  std::stringstream ss;
  std::string exported_name;
  if (parent_name != ""){
    exported_name = parent_name+ "_surface_classes";
  }
  else {
    exported_name = "surface_classes";
  }

  ss << exported_name << " = [\n";
  for (size_t i = 0; i < surface_classes.size(); i++) {
    const auto& item = surface_classes[i];
    if (i == 0) {
      ss << "    ";
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

std::string GenSubsystem::export_vec_elementary_molecule_types(std::ostream& out, PythonExportContext& ctx, const std::string& parent_name) {
  // prints vector into 'out' and returns name of the generated list
  std::stringstream ss;
  std::string exported_name;
  if (parent_name != ""){
    exported_name = parent_name+ "_elementary_molecule_types";
  }
  else {
    exported_name = "elementary_molecule_types";
  }

  ss << exported_name << " = [\n";
  for (size_t i = 0; i < elementary_molecule_types.size(); i++) {
    const auto& item = elementary_molecule_types[i];
    if (i == 0) {
      ss << "    ";
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

