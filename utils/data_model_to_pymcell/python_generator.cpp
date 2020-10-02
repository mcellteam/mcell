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

#include "generator_utils.h"
#include "python_generator.h"

using namespace std;
using namespace MCell::API;

namespace MCell {

using Json::Value;

void PythonGenerator::generate_single_parameter(std::ostream& out, Json::Value& parameter) {
  out << "# " << parameter[KEY_PAR_DESCRIPTION].asString() << "\n";
  out << parameter[KEY_PAR_NAME].asString() << " = " << parameter[KEY_PAR_EXPRESSION].asString();
  string units = parameter[KEY_PAR_UNITS].asString();
  if (units != "") {
    out << " # units: " << units;
  }
  out << "\n\n";
}


void PythonGenerator::generate_parameters(std::ostream& out) {
  Value& parameter_system = get_node(mcell, KEY_PARAMETER_SYSTEM);
  if (parameter_system.isMember(KEY_MODEL_PARAMETERS)) {
    Value& parameter_list = get_node(parameter_system, KEY_MODEL_PARAMETERS);
    for (Value::ArrayIndex i = 0; i < parameter_list.size(); i++) {
      generate_single_parameter(out, parameter_list[i]);
    }
  }
  out << "\n";
}


std::string PythonGenerator::generate_component_type(
    std::ostream& out, Json::Value& bngl_component_item, const std::string& mol_type_name) {

  string name = mol_type_name + '_' + make_id(bngl_component_item[KEY_CNAME].asString());
  gen_ctor_call(out, name, NAME_CLASS_COMPONENT_TYPE);

  vector<string> state_names;
  if (bngl_component_item.isMember(KEY_CSTATES)) {
    Value& cstates = bngl_component_item[KEY_CSTATES];
    for (Value::ArrayIndex i = 0; i < cstates.size(); i++) {
      state_names.push_back(cstates[i].asString());
    }
  }
  gen_param_list(out, NAME_STATES, state_names, false, true);
  out << CTOR_END;

  return name;
}


std::string PythonGenerator::generate_single_species_or_mol_type(
    std::ostream& out, Json::Value& molecule_list_item,
    const bool generate_species, const std::vector<string>& component_names) {

  assert(generate_species == component_names.empty());

  string name = make_id(molecule_list_item[KEY_MOL_NAME].asString());

  gen_ctor_call(out, name, (generate_species) ? NAME_CLASS_SPECIES : NAME_CLASS_ELEMENTARY_MOLECULE_TYPE);

  gen_param(out, NAME_NAME, molecule_list_item[KEY_MOL_NAME].asString(), true); // using original name

  if (!generate_species) {
    gen_param_list(out, NAME_COMPONENTS, component_names, true);
  }

  bool has_target_only = molecule_list_item[KEY_TARGET_ONLY].asBool();
  bool has_custom_time_step = molecule_list_item[KEY_CUSTOM_TIME_STEP].asString() != "";
  bool has_custom_space_step = molecule_list_item[KEY_CUSTOM_SPACE_STEP].asString() != "";
  bool has_extra_args = has_target_only || has_custom_time_step || has_custom_space_step;

  string mol_type = molecule_list_item[KEY_MOL_TYPE].asString();
  CHECK_PROPERTY(mol_type == VALUE_MOL_TYPE_2D || mol_type == VALUE_MOL_TYPE_3D);
  if (mol_type == VALUE_MOL_TYPE_3D) {
    gen_param_expr(out, NAME_DIFFUSION_CONSTANT_3D, molecule_list_item[KEY_DIFFUSION_CONSTANT], has_extra_args);
  }
  else {
    gen_param_expr(out, NAME_DIFFUSION_CONSTANT_2D, molecule_list_item[KEY_DIFFUSION_CONSTANT], has_extra_args);
  }

  release_assert(!(has_custom_time_step && has_custom_space_step) && "Only one of custom time or space step may be set");
  if (has_custom_time_step) {
    gen_param_expr(out, NAME_CUSTOM_TIME_STEP, molecule_list_item[KEY_CUSTOM_TIME_STEP], has_target_only);
  }
  else if (has_custom_space_step) {
    gen_param_expr(out, NAME_CUSTOM_SPACE_STEP, molecule_list_item[KEY_CUSTOM_SPACE_STEP], has_target_only);
  }

  if (has_target_only) {
    gen_param(out, NAME_TARGET_ONLY, has_target_only, false);
  }

  out << CTOR_END;

  return name;
}


SpeciesOrMolType PythonGenerator::generate_single_species_or_mol_type_w_components(
    std::ostream& out, Json::Value& molecule_list_item) {

  bool has_components = false;
  if (molecule_list_item.isMember(KEY_BNGL_COMPONENT_LIST) && molecule_list_item[KEY_BNGL_COMPONENT_LIST].size() > 0) {
    has_components = true;
  }

  // molecules in CellBlender represent either Species if they are a simple ComplexInstance or ElementaryMoleculeTypes

  // two different generation modes - when components are present then
  // we generate ElementaryMoleculeType, Species are a specific instantiation,
  // otherwise Species without components can be defined directly in the MCell3 species style
  if (has_components) {
    vector<string> component_names;

    string component_prefix = make_id(molecule_list_item[KEY_MOL_NAME].asString());

    // Components
    Value& bngl_component_list = get_node(molecule_list_item, KEY_BNGL_COMPONENT_LIST);
    for (Value::ArrayIndex i = 0; i < bngl_component_list.size(); i++) {
      string component_name = generate_component_type(out, bngl_component_list[i], component_prefix);
      component_names.push_back(component_name);
    }

    // Molecule Type
    string name = generate_single_species_or_mol_type(out, molecule_list_item, false, component_names);
    assert(name == component_prefix);

    return SpeciesOrMolType(name, false);
  }
  else {
    string name = generate_single_species_or_mol_type(out, molecule_list_item, true);

    return SpeciesOrMolType(name, true);
  }
}


void PythonGenerator::generate_species_and_mol_types(
    std::ostream& out, std::vector<SpeciesOrMolType>& species_and_mt_info) {

  Value& define_molecules = get_node(mcell, KEY_DEFINE_MOLECULES);
  check_version(KEY_DEFINE_MOLECULES, define_molecules, VER_DM_2014_10_24_1638);

  Value& molecule_list = get_node(define_molecules, KEY_MOLECULE_LIST);
  for (Value::ArrayIndex i = 0; i < molecule_list.size(); i++) {
    Value& molecule_list_item = molecule_list[i];
    check_version(KEY_MOLECULE_LIST, molecule_list_item, VER_DM_2018_10_16_1632);

    SpeciesOrMolType info = generate_single_species_or_mol_type_w_components(out, molecule_list_item);
    species_and_mt_info.push_back(info);
  }
}

} /* namespace MCell */
