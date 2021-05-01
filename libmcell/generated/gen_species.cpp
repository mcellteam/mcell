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
#include "gen_species.h"
#include "api/species.h"
#include "api/complex.h"
#include "api/elementary_molecule.h"
#include "api/species.h"

namespace MCell {
namespace API {

void GenSpecies::check_semantics() const {
}

void GenSpecies::set_initialized() {
  vec_set_initialized(elementary_molecules);
  initialized = true;
}

void GenSpecies::set_all_attributes_as_default_or_unset() {
  class_name = "Species";
  name = STR_UNSET;
  diffusion_constant_2d = FLT_UNSET;
  diffusion_constant_3d = FLT_UNSET;
  custom_time_step = FLT_UNSET;
  custom_space_step = FLT_UNSET;
  target_only = false;
  name = STR_UNSET;
  elementary_molecules = std::vector<std::shared_ptr<ElementaryMolecule>>();
  orientation = Orientation::DEFAULT;
  compartment_name = STR_UNSET;
}

bool GenSpecies::__eq__(const Species& other) const {
  return
    name == other.name &&
    diffusion_constant_2d == other.diffusion_constant_2d &&
    diffusion_constant_3d == other.diffusion_constant_3d &&
    custom_time_step == other.custom_time_step &&
    custom_space_step == other.custom_space_step &&
    target_only == other.target_only &&
    name == other.name &&
    vec_ptr_eq(elementary_molecules, other.elementary_molecules) &&
    orientation == other.orientation &&
    compartment_name == other.compartment_name;
}

bool GenSpecies::eq_nonarray_attributes(const Species& other, const bool ignore_name) const {
  return
    (ignore_name || name == other.name) &&
    diffusion_constant_2d == other.diffusion_constant_2d &&
    diffusion_constant_3d == other.diffusion_constant_3d &&
    custom_time_step == other.custom_time_step &&
    custom_space_step == other.custom_space_step &&
    target_only == other.target_only &&
    (ignore_name || name == other.name) &&
    true /*elementary_molecules*/ &&
    orientation == other.orientation &&
    compartment_name == other.compartment_name;
}

std::string GenSpecies::to_str(const std::string ind) const {
  std::stringstream ss;
  ss << get_object_name() << ": " <<
      "name=" << name << ", " <<
      "diffusion_constant_2d=" << diffusion_constant_2d << ", " <<
      "diffusion_constant_3d=" << diffusion_constant_3d << ", " <<
      "custom_time_step=" << custom_time_step << ", " <<
      "custom_space_step=" << custom_space_step << ", " <<
      "target_only=" << target_only << ", " <<
      "name=" << name << ", " <<
      "\n" << ind + "  " << "elementary_molecules=" << vec_ptr_to_str(elementary_molecules, ind + "  ") << ", " << "\n" << ind + "  " <<
      "orientation=" << orientation << ", " <<
      "compartment_name=" << compartment_name;
  return ss.str();
}

py::class_<Species> define_pybinding_Species(py::module& m) {
  return py::class_<Species, Complex, std::shared_ptr<Species>>(m, "Species", "There are three ways how to use this class:\n1) definition of simple species - in this case 'name' is \na single identifier and at least 'diffusion_constant_2d' or \n'diffusion_constant_3d' must be provided.\nExample: m.Species('A', diffusion_constant_3d=1e-6). \nSuch a definition must be added to subsystem or model so that  \nduring model initialization this species is transformed to MCell \nrepresentation and an ElementaryMoleculeType 'A' with a given \ndiffusion constant is created as well.\n2) full definition of complex species - in this case the \ninherited attribute 'elementary_molecules' from Complex\nis used as a definition of the complex and this gives information \non diffusion constants of the used elementary molecules.\nExample: m.Species(elementary_molecules=[ei1, ei2]). \nSuch a definition must be added to subsystem or model.   \n3) declaration of species - in this case only 'name' in the form of \nan BNGL string is provided. The complex instance specified by the name \nmust be fully qualified (i.e. all components are present and those \ncomponents that have a state have their state set).\nNo information on diffusion constants and other properties of \nused elementary molecules is provided, it must be provided elsewhere.\nExample: m.Species('A(b!1).B(a!1)').\nThis is a common form of usage when reaction rules are provided in a BNGL file.\nSuch declaration does no need to be added to subsystem or model.\nThis form is used as argument in cases where a fully qualified instance  \nmust be provided such as in molecule releases.\n \n")
      .def(
          py::init<
            const std::string&,
            const double,
            const double,
            const double,
            const double,
            const bool,
            const std::vector<std::shared_ptr<ElementaryMolecule>>,
            const Orientation,
            const std::string&
          >(),
          py::arg("name") = STR_UNSET,
          py::arg("diffusion_constant_2d") = FLT_UNSET,
          py::arg("diffusion_constant_3d") = FLT_UNSET,
          py::arg("custom_time_step") = FLT_UNSET,
          py::arg("custom_space_step") = FLT_UNSET,
          py::arg("target_only") = false,
          py::arg("elementary_molecules") = std::vector<std::shared_ptr<ElementaryMolecule>>(),
          py::arg("orientation") = Orientation::DEFAULT,
          py::arg("compartment_name") = STR_UNSET
      )
      .def("check_semantics", &Species::check_semantics)
      .def("__str__", &Species::to_str, py::arg("ind") = std::string(""))
      .def("__eq__", &Species::__eq__, py::arg("other"))
      .def("inst", &Species::inst, py::arg("orientation") = Orientation::DEFAULT, py::arg("compartment_name") = STR_UNSET, "Creates a copy of a Complex from this Species with specified orientation and compartment name. \n\n- orientation: Maximum one of orientation or compartment_name can be set, not both.\n\n- compartment_name: Maximum one of orientation or compartment_name can be set, not both.\n\n")
      .def("dump", &Species::dump)
      .def_property("name", &Species::get_name, &Species::set_name, "Name of the species in the BNGL format. \nOne must either specify name or elementary_molecules (inherited from Complex). \nThis argument name is parsed during model initialization.    \n")
      .def_property("diffusion_constant_2d", &Species::get_diffusion_constant_2d, &Species::set_diffusion_constant_2d, "This molecule is constrained to surface  with diffusion constant D. \nD can be zero, in which case the molecule doesn’t move. \nThe units of D are cm^2/s.\n")
      .def_property("diffusion_constant_3d", &Species::get_diffusion_constant_3d, &Species::set_diffusion_constant_3d, "This molecule diffuses in space with diffusion constant D. \nD can be zero, in which case the molecule doesn’t move. \nThe units of D are cm^2/s.\n \n")
      .def_property("custom_time_step", &Species::get_custom_time_step, &Species::set_custom_time_step, "Optional setting of a custom time step for this specific species. \nA molecule of this species should take timesteps of length custom_time_step (in seconds). \nUse either this or custom_time_step.\n")
      .def_property("custom_space_step", &Species::get_custom_space_step, &Species::set_custom_space_step, "Optional setting of a custom space step for this specific species. \nA molecule of this species should take steps of average length custom_space_step (in microns). \nUse either this or custom_time_step.\n     \n")
      .def_property("target_only", &Species::get_target_only, &Species::set_target_only, "A molecule of this species will not initiate reactions when it runs into other molecules. This\nsetting can speed up simulations when applied to a molecule at high concentrations \nthat reacts with a molecule at low concentrations (it is more efficient for\nthe low-concentration molecule to trigger the reactions). This directive does\nnot affect unimolecular reactions.\n")
    ;
}

std::string GenSpecies::export_to_python(std::ostream& out, PythonExportContext& ctx) {
  if (!export_even_if_already_exported() && ctx.already_exported(this)) {
    return ctx.get_exported_name(this);
  }
  std::string exported_name = std::string("species") + "_" + (is_set(name) ? fix_id(name) : std::to_string(ctx.postinc_counter("species")));
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
  ss << "m.Species(" << nl;
  if (name != STR_UNSET) {
    ss << ind << "name = " << "'" << name << "'" << "," << nl;
  }
  if (elementary_molecules != std::vector<std::shared_ptr<ElementaryMolecule>>() && !skip_vectors_export()) {
    ss << ind << "elementary_molecules = " << export_vec_elementary_molecules(out, ctx, exported_name) << "," << nl;
  }
  if (orientation != Orientation::DEFAULT) {
    ss << ind << "orientation = " << orientation << "," << nl;
  }
  if (compartment_name != STR_UNSET) {
    ss << ind << "compartment_name = " << "'" << compartment_name << "'" << "," << nl;
  }
  if (diffusion_constant_2d != FLT_UNSET) {
    ss << ind << "diffusion_constant_2d = " << f_to_str(diffusion_constant_2d) << "," << nl;
  }
  if (diffusion_constant_3d != FLT_UNSET) {
    ss << ind << "diffusion_constant_3d = " << f_to_str(diffusion_constant_3d) << "," << nl;
  }
  if (custom_time_step != FLT_UNSET) {
    ss << ind << "custom_time_step = " << f_to_str(custom_time_step) << "," << nl;
  }
  if (custom_space_step != FLT_UNSET) {
    ss << ind << "custom_space_step = " << f_to_str(custom_space_step) << "," << nl;
  }
  if (target_only != false) {
    ss << ind << "target_only = " << target_only << "," << nl;
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

std::string GenSpecies::export_vec_elementary_molecules(std::ostream& out, PythonExportContext& ctx, const std::string& parent_name) {
  // does not print the array itself to 'out' and returns the whole list
  std::stringstream ss;
  ss << "[";
  for (size_t i = 0; i < elementary_molecules.size(); i++) {
    const auto& item = elementary_molecules[i];
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

