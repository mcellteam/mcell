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
#include "libmcell/api/api_utils.h"
#include "bng/bng_defines.h"
#include "bng/bngl_names.h"

using namespace std;

namespace MCell {

using Json::Value;
using namespace API;

void BNGLGenerator::generate_single_bngl_parameter(Value& parameter) {
  if (parameter[KEY_PAR_DESCRIPTION].asString() != "") {
    bng_out << IND << "# " << parameter[KEY_PAR_DESCRIPTION].asString() << "\n";
  }
  bng_out << IND << fix_param_id(parameter[KEY_PAR_NAME].asString()) << " " <<
      replace_function_calls_in_expr(parameter[KEY_PAR_EXPRESSION].asString(), false);
  string units = parameter[KEY_PAR_UNITS].asString();
  if (units != "") {
    bng_out << " # units: " << units;
  }
  bng_out << "\n\n";
}


void BNGLGenerator::generate_single_python_parameter(std::ostream& python_out, Value& parameter) {
  string name = fix_param_id(parameter[KEY_PAR_NAME].asString());

  // skip MCELL_REDEFINE_ params
  if (name.find(BNG::MCELL_REDEFINE_PREFIX) == 0) {
    return;
  }

  data.check_if_already_defined_and_add(name, NAME_PARAMETER);
  python_out << name << " = " << VAR_BNGL_PARAMS << "['" << name << "']\n";
}


// NOTE: this belongs rather to the python generator
void BNGLGenerator::generate_parameters(std::ostream& python_out) {

  python_out << "# load parameters from BNGL\n";
  python_out <<
      VAR_BNGL_PARAMS << " = m.bngl_utils." << NAME_LOAD_BNGL_PARAMETERS << "(" <<
      get_abs_path(bngl_filename) << ", " <<
      S(SHARED) + "." + PARAMETER_OVERRIDES << ")\n\n";

  // and generate BNGL parameters and also their Python representations
  bng_out << BNG::BEGIN_PARAMETERS << "\n";
  Value& parameter_system = get_node(data.mcell, KEY_PARAMETER_SYSTEM);
  if (parameter_system.isMember(KEY_MODEL_PARAMETERS)) {
    Value& parameter_list = get_node(parameter_system, KEY_MODEL_PARAMETERS);
    for (Value::ArrayIndex i = 0; i < parameter_list.size(); i++) {
      generate_single_bngl_parameter(parameter_list[i]);
      generate_single_python_parameter(python_out, parameter_list[i]);
    }
  }
  bng_out << BNG::END_PARAMETERS << "\n\n";
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

      if (i + 1 != bngl_component_list.size()) {
        bng_out << ",";
      }
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

  python_out << "\n";
}


void BNGLGenerator::generate_mol_types(std::ostream& python_out) {
  bng_out << BNG::BEGIN_MOLECULE_TYPES << "\n";

  python_out <<
      "# set additional information about species and molecule types that cannot be stored in BNGL,\n"
      "# elementary molecule types are already in the subsystem after they were loaded from BNGL\n"
      "def " << SET_BNGL_MOLECULE_TYPES_INFO << "(subsystem):\n";

  Value& define_molecules = get_node(data.mcell, KEY_DEFINE_MOLECULES);
  check_version(KEY_DEFINE_MOLECULES, define_molecules, VER_DM_2014_10_24_1638);

  Value& molecule_list = get_node(define_molecules, KEY_MOLECULE_LIST);
  if (molecule_list.empty()) {
    python_out << IND4 << "pass # no molecule types are defined\n";
  }

  for (Value::ArrayIndex i = 0; i < molecule_list.size(); i++) {
    Value& molecule_list_item = molecule_list[i];
    check_version(KEY_MOLECULE_LIST, molecule_list_item, VER_DM_2018_10_16_1632);

    generate_bngl_mol_type(molecule_list_item);
    generate_python_mol_type_info(python_out, molecule_list_item);
  }
  bng_out << BNG::END_MOLECULE_TYPES << "\n\n";
}


void BNGLGenerator::generate_single_compartment(Json::Value& model_object) {
  const string& name = model_object[KEY_NAME].asString();
  const string& membrane_name = model_object[KEY_MEMBRANE_NAME].asString();
  const string& parent_object = model_object[KEY_PARENT_OBJECT].asString();

  // 3d compartment first
  bng_out << IND <<
      name << " " <<
      "3" << " " <<
      // volume does not have a correct value yet
      "1" << " " <<
      // parent is the membrane, may be unset
      membrane_name << "\n";

  // 2d compartment next
  if (membrane_name != "") {
    bng_out << IND <<
        membrane_name << " " <<
        "2" << " " <<
        // surface volume is ignored
        "0" << " " <<
        // parent object may be unset
        parent_object << "\n";
  }
}


static void add_parent_compartments_recursively(
    SharedGenData& data, Value& model_object_list, Value& model_object) {

  // simply keep inserting until we reach the top compartment

  const std::string& vol_comp = model_object[KEY_NAME].asString();
  assert(vol_comp != "");
  data.used_compartments.insert(vol_comp);

  const std::string& surf_comp = model_object[KEY_MEMBRANE_NAME].asString();
  if (surf_comp != "") {
    data.used_compartments.insert(surf_comp);
  }

  const std::string& parent_comp = model_object[KEY_PARENT_OBJECT].asString();
  if (parent_comp != "") {
    // find corresponding parent
    for (Value::ArrayIndex i = 0; i < model_object_list.size(); i++) {
      Value& parent_object = model_object_list[i];
      if (parent_object[KEY_NAME].asString() == parent_comp) {
        add_parent_compartments_recursively(data, model_object_list, parent_object);
      }
    }
  }
}


void BNGLGenerator::generate_compartments() {
  Value& model_objects = get_node(data.mcell, KEY_MODEL_OBJECTS);
  check_version(KEY_MODEL_OBJECTS, model_objects, VER_DM_2018_01_11_1330);

  Value& model_object_list = get_node(model_objects, KEY_MODEL_OBJECT_LIST);

  bool generate_compartments = false;
  for (Value::ArrayIndex i = 0; i < model_object_list.size(); i++) {
    Value& model_object = model_object_list[i];

    // generate the compartments that we need for rxns and recursively their parents
    if (data.is_used_compartment(model_object)) {
      generate_compartments = true;

      // recursively add compartment parents to used compartments, we need to generate them as well
      // because their children reference them
      add_parent_compartments_recursively(data, model_object_list, model_object);
    }
  }
  // do not generate empty section
  if (!generate_compartments) {
    return;
  }

  bng_out << BNG::BEGIN_COMPARTMENTS << "\n";
  bng_out <<
      IND << "# - volumes of compartments do not have correct values \n" <<
      IND << "# - this file is loaded through Subsystem.load_bngl_molecule_types_and_reaction_rules and\n" <<
      IND << "#   so compartments here are declared only so that the BNGL file can be parsed\n" <<
      IND << "# - compartments themselves are defined in the Python code using attributes of the GeometryObject\n";

  for (Value::ArrayIndex i = 0; i < model_object_list.size(); i++) {
    Value& model_object = model_object_list[i];

    // generate only the compartments that we need for rxns
    if (data.is_used_compartment(model_object)) {
      generate_single_compartment(model_object);
    }
  }

  bng_out << BNG::END_COMPARTMENTS << "\n\n";
}


static void fix_dots_in_simple_substances(vector<string>& substances) {
  for (string& s: substances) {
    s = fix_dots_in_simple_species(s);
  }
}


static bool has_in_out_compartments(const vector<string>& substances) {
  bool res = false;
  for (const string& s: substances) {
    size_t pos_in = s.find(S("@") + BNG::COMPARTMENT_NAME_IN);
    size_t pos_out = s.find(S("@") + BNG::COMPARTMENT_NAME_OUT);
    if (pos_in != string::npos || pos_out != string::npos) {
      return true;
    }
  }
  return false;
}


static void check_that_only_allowed_orientations_are_set(
    const vector<string>& orientations, const bool has_in_out_compartments) {

  // if IN/OUT is not used, volume molecules must not have orientation and
  // surface molecules must be UP
  // NOTE: this check can be improed with detection of what type of complex it is surface or volume
  for (const string& s: orientations) {
    release_assert(
        ((s == "" || s == "'") ||
         (has_in_out_compartments && (s == "'" || s == ","))) &&
        "Orientation in BNGL is not supported, should have been checked before");
  }
}


std::string BNGLGenerator::generate_single_reaction_rule(Json::Value& reaction_list_item, const bool generate_name) {
  string rxn_type = reaction_list_item[KEY_RXN_TYPE].asString();
  CHECK_PROPERTY(rxn_type == VALUE_IRREVERSIBLE || rxn_type == VALUE_REVERSIBLE);
  bool is_reversible = rxn_type == VALUE_REVERSIBLE;

  string name = get_rxn_id(reaction_list_item, data.unnamed_rxn_counter);

  bng_out << IND;
  if (generate_name) {
    // printing out name all the time would make the BNGL file hard to read
    bng_out << name << ": ";
  }

  vector<string> reac_substances;
  vector<string> reac_orientations;

  // reactants
  parse_rxn_rule_side(reaction_list_item[KEY_REACTANTS], reac_substances, reac_orientations);
  fix_dots_in_simple_substances(reac_substances);
  bool has_in_out = has_in_out_compartments(reac_substances); // in/out must be used in reactants to be allowed
  check_that_only_allowed_orientations_are_set(reac_orientations, has_in_out);
  for (size_t i = 0; i < reac_substances.size(); i++) {
    bng_out << reac_substances[i];
    if (i != reac_substances.size() - 1) {
      bng_out << " + ";
    }
  }
  bng_out << " ";

  bng_out << ((is_reversible) ? "<->" : "->") << " ";

  vector<string> prod_substances;
  vector<string> prod_orientations;

  // products
  parse_rxn_rule_side(reaction_list_item[KEY_PRODUCTS], prod_substances, prod_orientations);
  fix_dots_in_simple_substances(prod_substances);

  check_that_only_allowed_orientations_are_set(prod_orientations, has_in_out);
  for (size_t i = 0; i < prod_substances.size(); i++) {
    bng_out << prod_substances[i];
    if (i != prod_substances.size() - 1) {
      bng_out << " + ";
    }
  }
  bng_out << " ";

  // rates
  bng_out << reaction_list_item[KEY_FWD_RATE].asString();
  if (is_reversible) {
    bng_out << ", " << reaction_list_item[KEY_BKWD_RATE].asString();
  }
  bng_out << "\n";

  return name;
}


// rather limited for now
bool BNGLGenerator::can_express_count_with_bngl(
  const bool single_term,
  const bool rxn_not_mol,
  const std::string& where_to_count,
  const std::string& orientation,
  const std::string& multiplier_str,
  const std::string& rxn_step
) const {

  if (!single_term) {
    return false;
  }
  if (rxn_not_mol) {
    return false;
  }
  if (where_to_count != "" && where_to_count != VALUE_COUNT_LOCATION_WORLD) {
    return false;
  }
  if (orientation != "") {
    return false;
  }
  if (multiplier_str != "") {
    return false;
  }

  if (rxn_step != "") {
    // TODO: add test that uses parameters for both examined values
    double rxn_step_val;
    bool ok = get_parameter_value(data.mcell, rxn_step, rxn_step_val);
    if (!ok) {
      return false;
    }
    double time_step_val;
    const string& time_step = data.mcell[KEY_INITIALIZATION][KEY_TIME_STEP].asString();
    ok = get_parameter_value(data.mcell, time_step, time_step_val);
    if (!ok) {
      return false;
    }
    // time step and rxn step must be the same
    if (!cmp_eq(rxn_step_val, time_step_val)) {
      return false;
    }
  }

  return true;
}


void BNGLGenerator::generate_single_count(
    const std::string& observable_name,
    const std::string& what_to_count,
    const bool molecules_not_species) {

  bng_out << IND <<
      ((molecules_not_species) ? BNG::OBSERVABLE_MOLECULES : BNG::OBSERVABLE_SPECIES) << " " <<
      fix_id(observable_name) << " " <<
      what_to_count << "\n";
}


bool BNGLGenerator::can_express_release_with_bngl(Json::Value& release_site_item) {
  // for now only volume molecules, little limited,
  // mainly for data model generated by MCell4 reading a BNGL file

  // must be released into an object
  if (release_site_item[KEY_SHAPE].asString() != VALUE_OBJECT) {
    return false;
  }

  // must have a single compartment (probably too strict)
  bool multiple_compartments;
  string compartment = get_single_compartment(release_site_item[KEY_MOLECULE].asString(), &multiple_compartments);
  if (multiple_compartments) {
    return false;
  }

  // compartment must correspond to the object, or it must be empty and the object must the special default_comartment
  if (compartment == "" && release_site_item[KEY_OBJECT_EXPR].asString() != S(BNG::DEFAULT_COMPARTMENT_NAME) + "[ALL]") {
    return false;
  }

  if (compartment != "") {
    if (data.surface_to_volume_compartments_map.count(compartment) == 0) {
      if (release_site_item[KEY_OBJECT_EXPR].asString() != compartment + "[ALL]") {
        return false;
      }
    }
    else {
      // get volume compartment name if this is a surface compartment because that is what appears in the data model
      string vol_compartment = data.surface_to_volume_compartments_map[compartment];
      if (release_site_item[KEY_OBJECT_EXPR].asString() != vol_compartment + "[ALL]") {
        return false;
      }
    }
  }

  if (release_site_item[KEY_QUANTITY_TYPE].asString() != VALUE_NUMBER_TO_RELEASE) {
    return false;
  }
  if (release_site_item[KEY_RELEASE_PROBABILITY].asString() != "1") {
    return false;
  }
  if (release_site_item[KEY_STDDEV].asString() != "0") {
    return false;
  }
  if (release_site_item[KEY_ORIENT].asString() != ";" &&
      release_site_item[KEY_ORIENT].asString() != "'" ) {
    return false;
  }
  if (release_site_item[KEY_PATTERN].asString() != "") {
    return false;
  }

  return true;
}


void BNGLGenerator::generate_single_release_site(
    const std::string& bngl_cplx,
    const std::string& quantity) {
  bng_out << IND << bngl_cplx << " " << fix_param_id(quantity) <<  "\n";
}


} /* namespace MCell */
