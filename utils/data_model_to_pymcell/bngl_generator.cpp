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
#include "generator_structs.h"
#include "bngl_generator.h"

using namespace std;

namespace MCell {

using Json::Value;
using namespace API;

void BNGLGenerator::generate_single_bngl_parameter(Value& parameter) {
  bng_out << IND << "# " << parameter[KEY_PAR_DESCRIPTION].asString() << "\n";
  bng_out << IND << parameter[KEY_PAR_NAME].asString() << " " << parameter[KEY_PAR_EXPRESSION].asString();
  string units = parameter[KEY_PAR_UNITS].asString();
  if (units != "") {
    bng_out << " # units: " << units;
  }
  bng_out << "\n\n";
}


void BNGLGenerator::generate_single_python_parameter(std::ostream& python_out, Value& parameter) {
  string name = parameter[KEY_PAR_NAME].asString();
  python_out << name << " = " << VAR_BNGL_PARAMS << "['" << name << "']\n";
}


void BNGLGenerator::generate_parameters(std::ostream& python_out) {

  python_out << "# load parameters from BNGL\n";
  python_out << VAR_BNGL_PARAMS << " = m.bngl_utils.load_bngl_parameters('" << bngl_filename << "')\n\n";

  // and generate BNGL parameters and also their Python representations
  bng_out << "begin parameters\n";
  Value& parameter_system = get_node(data.mcell, KEY_PARAMETER_SYSTEM);
  if (parameter_system.isMember(KEY_MODEL_PARAMETERS)) {
    Value& parameter_list = get_node(parameter_system, KEY_MODEL_PARAMETERS);
    for (Value::ArrayIndex i = 0; i < parameter_list.size(); i++) {
      generate_single_bngl_parameter(parameter_list[i]);
      generate_single_python_parameter(python_out, parameter_list[i]);
    }
  }
  bng_out << "end parameters\n\n";
  python_out << "\n";
}


void BNGLGenerator::generate_bngl_mol_type(Json::Value& molecule_list_item) {

  string name = make_id(molecule_list_item[KEY_MOL_NAME].asString());

  bng_out << IND << name;

  bool has_components = false;
  if (molecule_list_item.isMember(KEY_BNGL_COMPONENT_LIST) && molecule_list_item[KEY_BNGL_COMPONENT_LIST].size() > 0) {
    has_components = true;
  }

  if (has_components) {
    bng_out << "(";
    // Components
    Value& bngl_component_list = get_node(molecule_list_item, KEY_BNGL_COMPONENT_LIST);
    for (Value::ArrayIndex i = 0; i < bngl_component_list.size(); i++) {
      Value& bngl_component = bngl_component_list[i];
      bng_out << bngl_component[KEY_CNAME].asString();

      Value& cstates = bngl_component[KEY_CSTATES];
      for (Value::ArrayIndex i = 0; i < cstates.size(); i++) {
        bng_out << "~" << cstates[i].asString();
      }

      gen_comma(bng_out, i, bngl_component_list);
    }
    bng_out << ")";
  }

  bng_out << "\n";
}


void BNGLGenerator::generate_python_mol_type_info(
    std::ostream& python_out, Json::Value& molecule_list_item) {

  string name = make_id(molecule_list_item[KEY_MOL_NAME].asString());

  python_out << IND4 <<
      name << " = subsystem." << NAME_FIND_ELEMENTARY_MOLECULE_TYPE << "('" << name << "')\n";
  python_out << IND4 << "assert " << name << ", \"Elementary molecule type '" + name + "' was not found\"\n";

  string mol_type = molecule_list_item[KEY_MOL_TYPE].asString();
  CHECK_PROPERTY(mol_type == VALUE_MOL_TYPE_2D || mol_type == VALUE_MOL_TYPE_3D);
  python_out << IND4;
  if (mol_type == VALUE_MOL_TYPE_3D) {
    gen_assign(python_out, name, NAME_DIFFUSION_CONSTANT_3D, molecule_list_item[KEY_DIFFUSION_CONSTANT].asString());
  }
  else {
    gen_assign(python_out, name, NAME_DIFFUSION_CONSTANT_2D, molecule_list_item[KEY_DIFFUSION_CONSTANT].asString());
  }

  bool has_custom_time_step = molecule_list_item[KEY_CUSTOM_TIME_STEP].asString() != "";
  bool has_custom_space_step = molecule_list_item[KEY_CUSTOM_SPACE_STEP].asString() != "";
  CHECK_PROPERTY(!(has_custom_time_step && has_custom_space_step) && "Only one of custom time or space step may be set");
  if (has_custom_time_step) {
    python_out << IND4;
    gen_assign(python_out, name, NAME_CUSTOM_TIME_STEP, molecule_list_item[KEY_CUSTOM_TIME_STEP].asString());
  }
  else if (has_custom_space_step) {
    python_out << IND4;
    gen_assign(python_out, name, NAME_CUSTOM_SPACE_STEP, molecule_list_item[NAME_CUSTOM_SPACE_STEP].asString());
  }

  if (molecule_list_item[KEY_TARGET_ONLY].asBool()) {
    python_out << IND4;
    gen_assign(python_out, name, NAME_TARGET_ONLY, true);
  }
}


void BNGLGenerator::generate_mol_types(std::ostream& python_out) {
  bng_out << "begin molecule types\n";

  python_out <<
      "# set additional information about species and molecule types that cannot be stored in BNGL,\n"
      "# elementary molecule types are already in the subsystem or model after they were loaded from BNGL\n"
      "def set_bngl_molecule_types_info(subsystem):\n";

  Value& define_molecules = get_node(data.mcell, KEY_DEFINE_MOLECULES);
  check_version(KEY_DEFINE_MOLECULES, define_molecules, VER_DM_2014_10_24_1638);

  Value& molecule_list = get_node(define_molecules, KEY_MOLECULE_LIST);
  for (Value::ArrayIndex i = 0; i < molecule_list.size(); i++) {
    Value& molecule_list_item = molecule_list[i];
    check_version(KEY_MOLECULE_LIST, molecule_list_item, VER_DM_2018_10_16_1632);

    generate_bngl_mol_type(molecule_list_item);
    generate_python_mol_type_info(python_out, molecule_list_item);
  }

  python_out << "\n";
  bng_out << "end molecule types\n\n";
}


std::string BNGLGenerator::generate_single_reaction_rule(Json::Value& reaction_list_item, const bool generate_name) {
  string rxn_type = reaction_list_item[KEY_RXN_TYPE].asString();
  CHECK_PROPERTY(rxn_type == VALUE_IRREVERSIBLE || rxn_type == VALUE_REVERSIBLE);
  bool is_reversible = rxn_type == VALUE_REVERSIBLE;

  // generate name only when needed
  string name = reaction_list_item[KEY_RXN_NAME].asString();
  if (name == "") {
    bool ok = convert_reaction_name(reaction_list_item[KEY_NAME].asString(), name);

    if (!ok) {
      name = UNNAMED_REACTION_RULE_PREFIX + to_string(data.unnamed_rxn_counter);
      data.unnamed_rxn_counter++;
    }
  }

  bng_out << IND;
  if (generate_name) {
    // printing out name all the time would make the BNGL file hard to read
    bng_out << name << ": ";
  }
  bng_out <<
      reaction_list_item[KEY_REACTANTS].asString() << " " <<
      ((is_reversible) ? "<->" : "->") << " " <<
      reaction_list_item[KEY_PRODUCTS].asString() << " " <<
      reaction_list_item[KEY_FWD_RATE].asString();

  if (is_reversible) {
    bng_out << ", " << reaction_list_item[KEY_BKWD_RATE].asString();
  }
  bng_out << "\n";

  return name;
}


void BNGLGenerator::generate_python_decl_bngl_rxn_rule(std::ostream& python_out, const std::string& name) {
  python_out << "# declaration of rxn rule defined in BNGL and used here in Python\n";
  python_out <<
    name << " == " << get_module_name_w_prefix(data.output_files_prefix, SUBSYSTEM) << "." <<
    NAME_FIND_REACTION_RULE << "('" << name << "')\n";
}


} /* namespace MCell */
