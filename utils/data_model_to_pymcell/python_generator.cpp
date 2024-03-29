/******************************************************************************
 *
 * Copyright (C) 2020 by
 * The Salk Institute for Biological Studies
 *
 * Use of this source code is governed by an MIT-style
 * license that can be found in the LICENSE file or at
 * https://opensource.org/licenses/MIT.
 *
******************************************************************************/

#include <fstream>
#include <regex>
#include <ctype.h>

#include "libmcell/api/api_utils.h"

#include "python_generator.h"
#include "generator_utils.h"
#include "mcell4_generator.h"
#include "bng/bng_defines.h"

using namespace std;
using namespace MCell::API;

namespace MCell {

using Json::Value;

void PythonGenerator::generate_single_parameter(std::ostream& out, Json::Value& parameter) {
  if (parameter[KEY_PAR_DESCRIPTION].asString() != "") {
    out << "# " << parameter[KEY_PAR_DESCRIPTION].asString() << "\n";
  }
  string python_expr;
  // replace operator ^ with operator **
  python_expr = regex_replace(parameter[KEY_PAR_EXPRESSION].asString(), regex("\\^"), "**");
  python_expr = replace_function_calls_in_expr(python_expr, true);
  string name = fix_param_id(parameter[KEY_PAR_NAME].asString());
  data.check_if_already_defined_and_add(name, NAME_PARAMETER);

  string ind;
  if (!data.not_overridable_python_params) {
    // allow params to be overridden
    out << "if " << NOT_DEFINED << "('" << name << "'):\n";
    ind = IND4;
  }
  out << ind << name << " = " << python_expr;
  string units = parameter[KEY_PAR_UNITS].asString();
  if (units != "") {
    out << " # units: " << units;
  }
  out << "\n\n";
}


static void get_used_ids(const string& expr, vector<string>& used_ids) {
  used_ids.clear();

  // same automaton but with different actions is used in replace_function_calls_in_expr
  enum state_t {
    START,
    IN_ID,
    IN_NUM
  };
  state_t state = START;
  string curr_id;
  for (size_t i = 0; i < expr.size(); i++) {
    char c = expr[i];
    switch (state) {
      case START:
        if (isalpha(c) || c == '_') {
          state = IN_ID;
          curr_id += c;
        }
        else if (isdigit(c)) {
          state = IN_NUM;
        }
        break;
      case IN_ID:
        if (isalnum(c) || c == '_') {
          curr_id += c;
        }
        else {
          if (mdl_functions_to_py_bngl_map.count(curr_id) == 0) {
            // count only IDs
            used_ids.push_back(curr_id);
          }
          curr_id = "";
          state = START;
        }
        break;
      case IN_NUM:
        if (isdigit(c) || c == '.' || c == 'e' || c == '+' || c == '-') {
          // ok
        }
        else {
          if (isalpha(c) || c == '_') {
            state = IN_ID;
            curr_id += c;
          }
          else {
            state = START;
          }
        }
        break;
    }
  }
  if (state == IN_ID) {
    if (mdl_functions_to_py_bngl_map.count(curr_id) == 0) {
      // count only IDs
      used_ids.push_back(curr_id);
    }
  }
}


struct ExprDepNode {
  size_t index;
  vector<ExprDepNode*> parents;
  vector<ExprDepNode*> children;
};


// simple class used mainly to delete all created nodes
class ExprDepNodeContainer {
public:
  ~ExprDepNodeContainer() {
    for (ExprDepNode* n: nodes) {
      delete n;
    }
    nodes.clear();
  }

  ExprDepNode* create_node(const size_t index) {
    ExprDepNode* res = new ExprDepNode();
    nodes.push_back(res);
    res->index = index;
    return res;
  }

  vector<ExprDepNode*> nodes;
};


// collects
static void traverse_given_level(
    const ExprDepNode* n, const uint level, vector<size_t>& level_ordering) {

  if (level == 0) {
    level_ordering.push_back(n->index);
  }
  else {
    // recursion stops if there are no children
    for (ExprDepNode* child: n->children) {
      traverse_given_level(child, level-1, level_ordering);
    }
  }
}


static void level_order_traversal(
    ExprDepNode* root, vector<size_t>& ordering) {

  vector<size_t> level_ordering;

  uint level = 1; // we are not counting the root node
  do {
    level_ordering.clear();
    traverse_given_level(root, level, level_ordering);

    ordering.insert(ordering.end(), level_ordering.begin(), level_ordering.end());

    level++;
  } while (!level_ordering.empty());
}


// not very efficient but the parameter systems in data models are not so huge
static void order_by_dependencies(const vector<ExprDepNode*>& nodes, vector<size_t>& ordering) {
  vector<bool> already_handled(nodes.size(), false);

  while (ordering.size() != nodes.size()) {

    for (size_t i = 0; i < nodes.size(); i++) {
      ExprDepNode* n = nodes[i];

      if (already_handled[n->index]) {
        continue;
      }

      bool all_parents_are_handled = true;
      for (auto& p: n->parents) {
        if (!already_handled[p->index]) {
          // skip, not ready for this parameter yet
          all_parents_are_handled = false;
          break;
        }
      }

      if (!all_parents_are_handled) {
        continue;
      }

      ordering.push_back(n->index);
      already_handled[n->index] = true;
    }
  }
}


static void define_parameter_ordering(Value& parameter_list, vector<size_t>& ordering) {
  ExprDepNodeContainer nodes;
  map<string, ExprDepNode*> defines;

  // first create all nodes
  for (Value::ArrayIndex i = 0; i < parameter_list.size(); i++) {
    const string& name = parameter_list[i][KEY_PAR_NAME].asString();
    defines[name] = nodes.create_node(i);
  }

  // then define dependencies
  for (Value::ArrayIndex i = 0; i < parameter_list.size(); i++) {
    const string& name = parameter_list[i][KEY_PAR_NAME].asString();
    const string& expr = parameter_list[i][KEY_PAR_EXPRESSION].asString();

    ExprDepNode* n = defines[name];
    vector<string> used_ids;

    // e.g. our parameter is a = b + c
    // parents of a are b and c (since
    // and child of b, resp. c, is a
    get_used_ids(expr, used_ids);
    for (const string& id: used_ids) {
      ExprDepNode* used_node = defines[id];
      // check function name
      if (used_node == nullptr) {
        assert(mdl_functions_to_py_bngl_map.count(id) == 0);
        ERROR("Did not find definition of identifier '" + id + "' when converting parameters, some MDL "
            "constants or functions are not supported."
        );
      }
      // make a bidirectional link
      n->parents.push_back(used_node);
      used_node->children.push_back(n);
    }
  }

  // now we might have multiple roots (nodes without uses), create just one
  /*ExprDepNode* root = nodes.create_node(INDEX_INVALID);
  for (auto node_it: defines) {
    if (node_it.second->parents.empty()) {
      root->children.push_back(node_it.second);
      node_it.second->parents.push_back(root);
    }
  }*/

  // order nodes by dependencies
  order_by_dependencies(nodes.nodes, ordering);

  // finally do a level-order traversal of the created dependency graph
  //level_order_traversal(root, ordering);
}


void PythonGenerator::generate_parameters(std::ostream& out) {
  Value& parameter_system = get_node(mcell, KEY_PARAMETER_SYSTEM);
  if (parameter_system.isMember(KEY_MODEL_PARAMETERS)) {
    Value& parameter_list = get_node(parameter_system, KEY_MODEL_PARAMETERS);

    vector<size_t> ordering;
    // sort parameters by their dependence
    define_parameter_ordering(parameter_list, ordering);

    // define_parameter_ordering can insert several parameters multiple times into ordering
    // TODO: make a better fix if needed, skipping duplicates for now
    set<size_t> already_generated;
    for (Value::ArrayIndex i: ordering) {
      if (already_generated.count(i) == 0) {
        generate_single_parameter(out, parameter_list[i]);
      }
      already_generated.insert(i);
    }
  }
  out << "\n";
}


std::string PythonGenerator::generate_component_type(
    std::ostream& out, Json::Value& bngl_component_item, const std::string& mol_type_name) {

  const string& cname = bngl_component_item[KEY_CNAME].asString();

  bool has_states =
      bngl_component_item.isMember(KEY_CSTATES) &&
      !bngl_component_item[KEY_CSTATES].empty();

  string name = COMPONENT_PREFIX + mol_type_name + '_' + make_id(cname);

  // generate only once
  // TODO: do we need to check here that all the declarations are the same?
  if (data.is_already_defined(name, NAME_CLASS_COMPONENT_TYPE)) {
    return name;
  }

  data.check_if_already_defined_and_add(name, NAME_CLASS_COMPONENT_TYPE);
  gen_ctor_call(out, name, NAME_CLASS_COMPONENT_TYPE);
  gen_param(out, NAME_NAME, cname, has_states);

  vector<string> state_names;
  if (has_states) {
    Value& cstates = bngl_component_item[KEY_CSTATES];
    for (Value::ArrayIndex i = 0; i < cstates.size(); i++) {
      state_names.push_back(cstates[i].asString());
    }
    gen_param_list(out, NAME_STATES, state_names, false, true);
  }
  out << CTOR_END;

  return name;
}


std::string PythonGenerator::generate_single_species_or_mol_type(
    std::ostream& out, Json::Value& molecule_list_item,
    const bool generate_species, const std::vector<string>& component_names) {

  assert(generate_species == component_names.empty());

  string orig_name = molecule_list_item[KEY_MOL_NAME].asString();
  string name = make_id(orig_name);

  const char* type = (generate_species) ? NAME_CLASS_SPECIES : NAME_CLASS_ELEMENTARY_MOLECULE_TYPE;

  data.check_if_already_defined_and_add(name, type);

  gen_description(out, molecule_list_item);

  string python_name = (generate_species) ? name : ELEMENTARY_MOLECULE_TYPE_PREFIX + name;
  gen_ctor_call(out, python_name, type);

  string name_to_generate = orig_name;
  if (API::is_simple_species(orig_name)) {
    // simple species must not contain dots in their name
    name_to_generate = name;
  }

  gen_param(out, NAME_NAME, name_to_generate, true); // using original name

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

  CHECK_PROPERTY(!(has_custom_time_step && has_custom_space_step) && "Only one of custom time or space step may be set");
  if (has_custom_time_step) {
    gen_param_expr(out, NAME_CUSTOM_TIME_STEP, molecule_list_item[KEY_CUSTOM_TIME_STEP], has_target_only);
  }
  else if (has_custom_space_step) {
    gen_param_expr(out, NAME_CUSTOM_SPACE_STEP, molecule_list_item[KEY_CUSTOM_SPACE_STEP], has_target_only);
  }

  if (has_target_only) {
    gen_param(out, NAME_TARGET_ONLY, true, false);
  }

  out << CTOR_END;

  return python_name;
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


void PythonGenerator::get_surface_class_property_info(
    const string& sc_name,
    Value& property,
    string& name, string& type_name,
    string& affected_mols, string& orientation, string& clamp_concentration
) {
  check_version(KEY_SURFACE_CLASS_PROP_LIST, property, VER_DM_2015_11_08_1756);

  affected_mols = property[KEY_AFFECTED_MOLS].asString();
  if (affected_mols == BNG::ALL_MOLECULES) {
    affected_mols = S(MDOT) + NAME_CV_AllMolecules;
  }
  else if (affected_mols == BNG::ALL_VOLUME_MOLECULES) {
    affected_mols = S(MDOT) + NAME_CV_AllVolumeMolecules;
  }
  else if (affected_mols == BNG::ALL_SURFACE_MOLECULES) {
    affected_mols = S(MDOT) + NAME_CV_AllSurfaceMolecules;
  }
  else if (affected_mols == VALUE_SINGLE) {
    affected_mols = property[KEY_MOLECULE].asString();
  }

  orientation = convert_orientation(property[KEY_SURF_CLASS_ORIENT].asString());

  clamp_concentration = "";
  string surf_class_type = property[KEY_SURF_CLASS_TYPE].asString();
  if (surf_class_type == VALUE_TRANSPARENT) {
    type_name = NAME_EV_TRANSPARENT;
  }
  else if (surf_class_type == VALUE_REFLECTIVE) {
    type_name = NAME_EV_REFLECTIVE;
  }
  else if (surf_class_type == VALUE_ABSORPTIVE) {
    type_name = NAME_EV_ABSORPTIVE;
  }
  else if (surf_class_type == VALUE_CLAMP_CONCENTRATION) {
    type_name = NAME_EV_CONCENTRATION_CLAMP;
    clamp_concentration = property[KEY_CLAMP_VALUE].asString();
  }
  else if (surf_class_type == VALUE_CLAMP_FLUX) {
    type_name = NAME_EV_FLUX_CLAMP;
    clamp_concentration = property[KEY_CLAMP_VALUE].asString();
  }
  else if (surf_class_type == "") {
    type_name = NAME_EV_REACTIVE;
  }
  else {
    ERROR(S("Invalid type/property ") + KEY_SURF_CLASS_TYPE + " " + surf_class_type +
        " of surface class '" + sc_name + "0'.");
  }

  name = make_id(property[KEY_NAME].asString());
  if (name == "") {
    string tn_lower = type_name;
    transform(tn_lower.begin(), tn_lower.end(), tn_lower.begin(), ::tolower);
    // to be sure that there is no conflict, let's put a number at the end
    name = SURFACE_CLASS_PREFIX + tn_lower + "_" + affected_mols + to_string(unnamed_surf_class_counter);
    unnamed_surf_class_counter++;
  }
}


void PythonGenerator::generate_surface_classes(
    std::ostream& out, std::vector<std::string>& sc_names) {

  if (!mcell.isMember(KEY_DEFINE_SURFACE_CLASSES)) {
    return;
  }

  Value& define_surface_classes = get_node(mcell, KEY_DEFINE_SURFACE_CLASSES);
  check_version(KEY_DEFINE_SURFACE_CLASSES, define_surface_classes, VER_DM_2014_10_24_1638);
  Value& surface_class_list = get_node(define_surface_classes, KEY_SURFACE_CLASS_LIST);

  for (Value::ArrayIndex i = 0; i < surface_class_list.size(); i++) {
    Value& surface_class_list_item = surface_class_list[i];

    Value& surface_class_prop_list = get_node(surface_class_list_item, KEY_SURFACE_CLASS_PROP_LIST);

    string sc_name = make_id(surface_class_list_item[KEY_NAME].asString());
    sc_names.push_back(sc_name);

    vector<string> sc_prop_names;
    if (surface_class_prop_list.size() > 1) {
      for (Value::ArrayIndex i = 0; i < surface_class_prop_list.size(); i++) {
        Value& surface_class_prop_item = surface_class_prop_list[i];

        string name, type_name, affected_mols, orientation_name, clamp_concentration;
        get_surface_class_property_info(
            sc_name, surface_class_prop_item,
            name, type_name, affected_mols, orientation_name, clamp_concentration);

        bool is_clamp = type_name == NAME_EV_CONCENTRATION_CLAMP || type_name == NAME_EV_FLUX_CLAMP;
        bool has_affected_mols = affected_mols != "";

        sc_prop_names.push_back(name);
        data.check_if_already_defined_and_add(name, NAME_CLASS_SURFACE_PROPERTY);
        gen_ctor_call(out, name, NAME_CLASS_SURFACE_PROPERTY, true);
        gen_param_enum(out, NAME_TYPE, NAME_ENUM_SURFACE_PROPERTY_TYPE, type_name, has_affected_mols || is_clamp);

        if (has_affected_mols) {
          gen_param_expr(
              out, NAME_AFFECTED_COMPLEX_PATTERN,
              make_species_or_cplx(data, affected_mols, orientation_name), is_clamp);
        }

        if (is_clamp) {
          gen_param_expr(out, NAME_CONCENTRATION, clamp_concentration, false);
        }

        out << CTOR_END;
      }
    }

    gen_description(out, surface_class_list_item);

    data.check_if_already_defined_and_add(sc_name, NAME_CLASS_SURFACE_CLASS);
    gen_ctor_call(out, sc_name, NAME_CLASS_SURFACE_CLASS, true);
    gen_param(out, NAME_NAME, sc_name, true);

    if (!sc_prop_names.empty()) {
      // use a list of properties
      gen_param_list(out, NAME_PROPERTIES, sc_prop_names, false);
    }
    else {
      // simplified setup, directly set members
      string name, type_name, affected_mols, orientation_name, clamp_concentration;
      get_surface_class_property_info(
          sc_name, surface_class_prop_list[0],
          name, type_name, affected_mols, orientation_name, clamp_concentration);

      bool is_clamp = type_name == NAME_EV_CONCENTRATION_CLAMP || type_name == NAME_EV_FLUX_CLAMP;
      bool has_affected_mols = affected_mols != "";

      gen_param_enum(out, NAME_TYPE, NAME_ENUM_SURFACE_PROPERTY_TYPE, type_name, has_affected_mols || is_clamp);

      if (has_affected_mols) {
        gen_param_expr(
            out, NAME_AFFECTED_COMPLEX_PATTERN,
            make_species_or_cplx(data, affected_mols, orientation_name), is_clamp);
      }

      if (is_clamp) {
        gen_param_expr(out, NAME_CONCENTRATION, clamp_concentration, false);
      }
    }
    out << CTOR_END;
  }
}


void PythonGenerator::generate_variable_rate(const std::string& rate_array_name, Json::Value& variable_rate_text) {
  // append to the parameters file
  ofstream out;
  open_and_check_file_w_prefix(data.output_files_prefix, PARAMETERS, out, true);

  out << "\n" << make_section_comment("variable rate");
  out << rate_array_name << " = [\n";

  string vr = variable_rate_text.asString();
  size_t pos = 0;
  bool print_comma = false;
  while (pos < vr.size()) {

    size_t tab = vr.find('\t', pos + 1);
    size_t space = vr.find(' ', pos + 1);
    if (tab == string::npos || space < tab) {
      tab = space;
    }

    size_t nl = vr.find('\n', pos + 1);
    if (nl == string::npos) {
      nl = vr.size();
    }
    if (tab > nl || tab == pos || nl == pos) {
      ERROR("Malformed variable_rate_text in datamodel.");
    }
    if (tab == string::npos) {
      break; // end of input
    }

    string time = vr.substr(pos, tab - pos);
    string rate = vr.substr(tab + 1, nl - (tab + 1));
    time = trim(time);
    rate = trim(rate);

    if (print_comma) {
      out << ",\n";
    }
    print_comma = true;

    out << "  [" << time << ", " << rate << "]";

    pos = nl;
    pos++;
  }
  out << "\n] # " << rate_array_name << "\n\n";

  out.close();
}



void PythonGenerator::generate_rxn_rule_side(std::ostream& out, Json::Value& substances_node) {

  string str = substances_node.asString();

  // special case for rxns without products
  if (str == NULL_PRODUCTS) {
    out << "[ ]";
    return;
  }

  vector<string> substances;
  vector<string> orientations;
  parse_rxn_rule_side(substances_node, substances, orientations);

  out << "[ ";
  for (size_t i = 0; i < substances.size(); i++) {

    string orient = convert_orientation(orientations[i], true);
    string compartment = get_single_compartment(substances[i]);
    out << make_species_or_cplx(data, substances[i], orient, compartment);
    print_comma(out, i, substances);
  }
  out << " ]";
}


std::string PythonGenerator::generate_single_reaction_rule(std::ostream& out, Json::Value& reaction_list_item) {
  check_version(KEY_MOLECULE_LIST, reaction_list_item, VER_DM_2018_01_11_1330);

  string name = get_rxn_id(reaction_list_item, data.unnamed_rxn_counter);

  gen_description(out, reaction_list_item);

  data.check_if_already_defined_and_add(name, NAME_CLASS_REACTION_RULE);
  gen_ctor_call(out, name, NAME_CLASS_REACTION_RULE);
  gen_param(out, NAME_NAME, name, true);

  // single line for now
  out << IND << NAME_REACTANTS << " = ";
  generate_rxn_rule_side(out, reaction_list_item[KEY_REACTANTS]);
  out << ",\n";

  out << IND << NAME_PRODUCTS << " = ";
  generate_rxn_rule_side(out, reaction_list_item[KEY_PRODUCTS]);
  out << ",\n";

  if (reaction_list_item[KEY_VARIABLE_RATE_SWITCH].asBool()) {
    // variable rates
    CHECK_PROPERTY(reaction_list_item[KEY_VARIABLE_RATE_VALID].asBool() && "variable_rate_switch must be equal to variable_rate_valid");
    string rate_array_name = reaction_list_item[KEY_VARIABLE_RATE].asString();
    size_t dot = rate_array_name.rfind('.');
    if (dot != string::npos) {
      // remove file extension
      rate_array_name = rate_array_name.substr(0, dot);
    }
    // and remove possibly other dots in the name
    rate_array_name = make_id(rate_array_name);
    generate_variable_rate(rate_array_name, reaction_list_item[KEY_VARIABLE_RATE_TEXT]);
    gen_param_id(out, NAME_VARIABLE_RATE, rate_array_name, false); // module parameters is imported as *
  }
  else {
    // fwd or rev rates
    string rxn_type = reaction_list_item[KEY_RXN_TYPE].asString();
    CHECK_PROPERTY(rxn_type == VALUE_IRREVERSIBLE || rxn_type == VALUE_REVERSIBLE);
    bool is_reversible = rxn_type == VALUE_REVERSIBLE;

    gen_param_expr(out, NAME_FWD_RATE, reaction_list_item[KEY_FWD_RATE], is_reversible);
    if (is_reversible) {
      gen_param(out, NAME_REV_NAME, name + REV_RXN_SUFFIX, true);
      gen_param_expr(out, NAME_REV_RATE, reaction_list_item[KEY_BKWD_RATE], false);
    }
  }

  out << CTOR_END;

  return name;
}


void PythonGenerator::generate_reaction_rules(
    ostream& out, std::vector<IdLoc>& rxn_names) {

  Value& define_reactions = get_node(mcell, KEY_DEFINE_REACTIONS);
  check_version(KEY_DEFINE_REACTIONS, define_reactions, VER_DM_2014_10_24_1638);

  Value& reaction_list = get_node(define_reactions, KEY_REACTION_LIST);
  for (Value::ArrayIndex i = 0; i < reaction_list.size(); i++) {
    Value& reaction_list_item = reaction_list[i];
    string name = generate_single_reaction_rule(out, reaction_list_item);
    rxn_names.push_back(IdLoc(name, true));
  }
}


string PythonGenerator::generate_single_geometry_object(
    ostream& out, const int index, Value& object) {

  string parent_name = S(KEY_OBJECT_LIST) + "[" + to_string(index) + "]";

  string name = make_id(get_node(parent_name, object, KEY_NAME).asString());
  data.check_if_already_defined_and_add(name, NAME_CLASS_GEOMETRY_OBJECT);

  Value& vertex_list = get_node(parent_name, object, KEY_VERTEX_LIST);
  Value& element_connections = get_node(parent_name, object, KEY_ELEMENT_CONNECTIONS);
  // TODO: material_names

  out << make_start_block_comment(name);

  // vertex_list
  string id_vertex_list = name + "_" + NAME_VERTEX_LIST;
  out << id_vertex_list << " = [\n";
  for (Value::ArrayIndex i = 0; i < vertex_list.size(); i++) {
    out << IND << "[";
    Value& vertex = vertex_list[i];
    for (Value::ArrayIndex k = 0; k < vertex.size(); k++) {
      out << vertex[k].asDouble();
      gen_comma(out, k, vertex);
    }
    out << "]";
    gen_comma(out, i, vertex_list);
    out << "\n";
  }
  out << "] # " << id_vertex_list << "\n\n";

  // element_connections
  string id_element_connections = name + "_" + NAME_WALL_LIST;
  out << id_element_connections << " = [\n";
  for (Value::ArrayIndex i = 0; i < element_connections.size(); i++) {
    out << IND << "[";
    Value& element = element_connections[i];
    for (Value::ArrayIndex k = 0; k < element.size(); k++) {
      out << element[k].asInt();
      gen_comma(out, k, element);
    }
    out << "]";
    gen_comma(out, i, element_connections);
    out << "\n";
  }
  out << "] # " << id_element_connections << "\n\n";

  // surface areas
  vector<string> sr_global_names;
  if (object.isMember(KEY_DEFINE_SURFACE_REGIONS)) {
    Value& define_surface_regions = object[KEY_DEFINE_SURFACE_REGIONS];
    for (Value::ArrayIndex i = 0; i < define_surface_regions.size(); i++) {

      string sr_name = make_id(get_node(
          KEY_DEFINE_SURFACE_REGIONS, define_surface_regions[i], KEY_NAME).asString());
      Value& include_elements = get_node(
          KEY_DEFINE_SURFACE_REGIONS, define_surface_regions[i], KEY_INCLUDE_ELEMENTS);

      string sr_global_name = name + "_" + sr_name;
      string sr_element_connections_name = sr_global_name + "_" + NAME_WALL_INDICES;
      out << sr_element_connections_name << " = [";
      for (Value::ArrayIndex k = 0; k < include_elements.size(); k++) {

        if (k % 16 == 0) {
          out << "\n" << IND;
        }
        out << include_elements[k].asInt();
        gen_comma(out, k, include_elements);
      }
      out << "\n] #" << sr_element_connections_name << "\n\n";

      out << sr_global_name << " = " << MDOT << NAME_CLASS_SURFACE_REGION << "(\n";
      out << IND << NAME_NAME << " = '" << sr_name << "',\n";
      out << IND << NAME_WALL_INDICES << " = " << sr_element_connections_name << "\n";
      out << CTOR_END;

      sr_global_names.push_back(sr_global_name);
    }
  }

  // object creation itself
  out << name << " = " << MDOT << NAME_CLASS_GEOMETRY_OBJECT << "(\n";
  gen_param(out, NAME_NAME, name, true);
  gen_param_id(out, NAME_VERTEX_LIST, id_vertex_list, true);
  gen_param_id(out, NAME_WALL_LIST, id_element_connections, true);
  out << IND << NAME_SURFACE_REGIONS << " = [";
  for (size_t i = 0; i < sr_global_names.size(); i++) {
    out << sr_global_names[i];
    print_comma(out, i, sr_global_names);
  }
  out << "]\n)\n";

  out << make_end_block_comment(name);

  return name;
}


void PythonGenerator::generate_geometry(std::ostream& out, std::vector<std::string>& geometry_objects) {

  Value& geometrical_objects = get_node(mcell, KEY_GEOMETRICAL_OBJECTS);
  if (!geometrical_objects.isMember(KEY_OBJECT_LIST)) {
    return;
  }
  Value& object_list = get_node(geometrical_objects, KEY_OBJECT_LIST);
  for (Value::ArrayIndex i = 0; i < object_list.size(); i++) {
    Value& object = object_list[i];
    string name = generate_single_geometry_object(out, i, object);
    if (name == BNG::DEFAULT_COMPARTMENT_NAME) {
      data.has_default_compartment_object = true;
    }
    geometry_objects.push_back(name);
  }
}


// returns true if the two release sites can be merged into a single one
static bool can_be_in_same_list_release_site(Value& rs1, Value& rs2) {
  return
      rs1[KEY_SHAPE].asString() == VALUE_LIST &&
      rs2[KEY_SHAPE].asString() == VALUE_LIST &&
      rs1[KEY_PATTERN].asString() == rs2[KEY_PATTERN].asString() &&
      rs1[KEY_QUANTITY].asString() == rs2[KEY_QUANTITY].asString() &&
      rs1[KEY_QUANTITY_TYPE].asString() == rs2[KEY_QUANTITY_TYPE].asString() &&
      rs1[KEY_RELEASE_PROBABILITY].asString() == rs2[KEY_RELEASE_PROBABILITY].asString() &&
      rs1[KEY_SITE_DIAMETER].asString() == rs2[KEY_SITE_DIAMETER].asString() &&
      rs1[KEY_STDDEV].asString() == rs2[KEY_STDDEV].asString();
}


std::string PythonGenerator::generate_single_molecule_release_info_array(
    std::ostream& out,
    std::string& rel_site_name,
    Json::Value& release_site_list,
    Json::Value::ArrayIndex begin,
    Json::Value::ArrayIndex end) {

  string name = MOLECULE_LIST_PREFIX + rel_site_name;

  out << name << " = [\n";
  for (Value::ArrayIndex rs_index = begin; rs_index < end; rs_index++) {
    Value& release_site_item = release_site_list[rs_index];

    Value& points_list = release_site_item[KEY_POINTS_LIST];

    for (Value::ArrayIndex i = 0; i < points_list.size(); i++) {
      out << "    " << MDOT << NAME_CLASS_MOLECULE_RELEASE_INFO << "(\n";

      string cplx = release_site_item[KEY_MOLECULE].asString();
      bool is_vol = is_volume_species(mcell, cplx);
      string orient = convert_orientation(release_site_item[KEY_ORIENT].asString(), !is_vol);
      out << "    ";
      gen_param_expr(out, NAME_COMPLEX,
          make_species_or_cplx(data, cplx, orient),
          true);

      Value& point = points_list[i];
      if (point.size() != 3) {
        ERROR("Release site " + rel_site_name + ": points_list item does not have three values.");
      }
      out << "        " << NAME_LOCATION << " = [" << point[0].asDouble() << ", " << point[1].asDouble() << ", " << point[2].asDouble() << "]";

      out << "\n    )";
      if (i + 1 != points_list.size()) {
        out << ", ";
      }
    }

    if (rs_index + 1 != end) {
      out << ", ";
    }
    out << "\n";
  }
  out << "] # " << name << "\n\n";

  return name;
}


static void gen_region_expr_assignment_for_rel_site(ostream& out, string region_expr) {
  // we don't really need to parse the expression, just dump it in a right way
  // example: Cell[ALL] - (Organelle_1[ALL] + Organelle_2[ALL])
  // remove [ALL]
  // replace [text] with _text

  regex pattern_all("\\[ALL\\]");
  region_expr = regex_replace(region_expr, pattern_all, "");

  regex pattern_surf_region("\\[([^\\]]*)\\]");
  region_expr = regex_replace(region_expr, pattern_surf_region, "_$1");

  gen_param_expr(out, NAME_REGION, region_expr, true);
}


static void error_release_pattern(const string& name) {
  ERROR("Requested release pattern " + name + " not found in data model.");
}


void PythonGenerator::generate_release_pattern(std::ostream& out, const std::string& name, std::string& delay_string) {

  // generate only once, however we must return the delay string
  bool already_generated = data.is_already_defined(name, NAME_CLASS_RELEASE_PATTERN);

  if (!mcell.isMember(KEY_DEFINE_RELEASE_PATTERNS)) {
    error_release_pattern(name);
  }
  Value& define_release_patterns = get_node(mcell, KEY_DEFINE_RELEASE_PATTERNS);
  Value& release_pattern_list = get_node(define_release_patterns, KEY_RELEASE_PATTERN_LIST);

  for (Value::ArrayIndex i = 0; i < release_pattern_list.size(); i++) {
    Value& release_pattern_item = release_pattern_list[i];
    if (release_pattern_item[KEY_NAME].asString() != name) {
      continue;
    }

    if (!already_generated) {
      // we found the release pattern
      check_version(KEY_RELEASE_PATTERN_LIST, release_pattern_item, VER_DM_2018_01_11_1330);

      data.check_if_already_defined_and_add(name, NAME_CLASS_RELEASE_PATTERN);
      gen_ctor_call(out, name, NAME_CLASS_RELEASE_PATTERN);
      gen_param(out, NAME_NAME, name, true);
      gen_param_expr(out, NAME_RELEASE_INTERVAL, release_pattern_item[KEY_RELEASE_INTERVAL], true);
      gen_param_expr(out, NAME_TRAIN_DURATION, release_pattern_item[KEY_TRAIN_DURATION], true);
      gen_param_expr(out, NAME_TRAIN_INTERVAL, release_pattern_item[KEY_TRAIN_INTERVAL], true);
      gen_param_expr(out, NAME_NUMBER_OF_TRAINS, release_pattern_item[KEY_NUMBER_OF_TRAINS], false);
      out << CTOR_END;
    }

    delay_string = release_pattern_item[KEY_DELAY].asString();
    return;
  }

  error_release_pattern(name);
}


void PythonGenerator::generate_release_sites(std::ostream& out, std::vector<std::string>& release_site_names) {

  if (!mcell.isMember(KEY_RELEASE_SITES)) {
    return;
  }

  Value& release_sites = mcell[KEY_RELEASE_SITES];
  check_version(KEY_RELEASE_SITES, release_sites, VER_DM_2014_10_24_1638);
  Value& release_site_list = release_sites[KEY_RELEASE_SITE_LIST];

  Value::ArrayIndex i = 0;
  while (i < release_site_list.size()) {
    Value& release_site_item = release_site_list[i];
    check_version(KEY_RELEASE_SITE_LIST, release_site_item, VER_DM_2018_01_11_1330);

    string name = make_id(release_site_item[KEY_NAME].asString());
    string shape = release_site_item[KEY_SHAPE].asString();
    string molecule_list_name = "";
    if (shape == VALUE_LIST) {
      // 1) try to find the largest number of subsequent
      // release sites that can be represented by a single ReleaseSite
      Value::ArrayIndex matching_end = i + 1;
      while (matching_end < release_site_list.size() &&
          can_be_in_same_list_release_site(release_site_item, release_site_list[matching_end])) {
        matching_end++;
      }

      // remove suffix from release name
      if (i < matching_end - 1 && name.substr(name.size() - 2) == "_0") {
        name = name.substr(0, name.size() - 2);
      }

      // 2) generate an array of SingleMoleculeReleaseInfo objects
      molecule_list_name =
          generate_single_molecule_release_info_array(out, name, release_site_list, i, matching_end);

      // skip all release sites we handled here
      i = matching_end - 1;
    }

    // generate release pattern if needed
    string rel_pat_name = release_site_item[KEY_PATTERN].asString();
    string delay_string = "";
    if (rel_pat_name != "") {
      generate_release_pattern(out, rel_pat_name, delay_string);
    }

    release_site_names.push_back(name);

    gen_description(out, release_site_item);

    data.check_if_already_defined_and_add(name, NAME_CLASS_RELEASE_SITE);
    gen_ctor_call(out, name, NAME_CLASS_RELEASE_SITE);
    gen_param(out, NAME_NAME, name, true);

    string compartment;
    if (shape != VALUE_LIST) {
      string cplx = release_site_item[KEY_MOLECULE].asString();
      if (cplx == "") {
        ERROR("Release site '" + name + "' does not have its molecule/complex to be released set.");
      }
      bool is_vol = is_volume_species(data.mcell, cplx);
      string orientation = convert_orientation(release_site_item[KEY_ORIENT].asString(), !is_vol);
      if (orientation != "" && is_vol) {
          cout <<
              "Ignoring orientation set for release site " << name << " with species " << cplx <<
              ", this species represent volume molecules.\n";
          orientation = "";
      }
      compartment = get_single_compartment(cplx);
      gen_param_expr(out, NAME_COMPLEX, make_species_or_cplx(data, cplx, orientation, compartment), true);
    }

    if (delay_string != "" && delay_string != "0") {
      gen_param_expr(out, NAME_RELEASE_TIME, delay_string, true);
    }

    if (rel_pat_name != "") {
      gen_param_id(out, NAME_RELEASE_PATTERN, rel_pat_name, true);
    }

    string prob = release_site_item[KEY_RELEASE_PROBABILITY].asString();
    bool prob_not_1 = (prob != "1" && prob != "1.0");

    if (shape == VALUE_SPHERICAL) {
      if (compartment != "") {
        ERROR("Cannot use compartment " + compartment + " with spherical release.");
      }
      gen_param_enum(out, NAME_SHAPE, NAME_ENUM_SHAPE, NAME_EV_SPHERICAL, true);
      gen_param_list_3_floats(
          out, NAME_LOCATION,
          release_site_item[KEY_LOCATION_X], release_site_item[KEY_LOCATION_Y], release_site_item[KEY_LOCATION_Z],
          true
      );
      gen_param_expr(out, NAME_SITE_DIAMETER, release_site_item[KEY_SITE_DIAMETER], true);
    }
    else if (shape == VALUE_OBJECT) {
      gen_region_expr_assignment_for_rel_site(out, release_site_item[KEY_OBJECT_EXPR].asString());
    }
    else if (shape == VALUE_LIST) {
      // FIXME: check that no compartments are used here
      assert(molecule_list_name != "");
      bool diam_is_zero = release_site_item[KEY_SITE_DIAMETER] == "0";
      gen_param_id(out, NAME_MOLECULE_LIST, molecule_list_name, !diam_is_zero);
      if (!diam_is_zero) {
        gen_param_expr(out, NAME_SITE_DIAMETER, release_site_item[KEY_SITE_DIAMETER], prob_not_1);
      }
    }
    else {
      ERROR("Shape " + shape + " is not supported yet");
    }

    if (shape != VALUE_LIST) {
      string quantity_type = release_site_item[KEY_QUANTITY_TYPE].asString();
      if (quantity_type == VALUE_NUMBER_TO_RELEASE) {
        check_not_empty(release_site_item, KEY_QUANTITY, "Release site");
        gen_param_expr(out, NAME_NUMBER_TO_RELEASE, release_site_item[KEY_QUANTITY], prob_not_1);
      }
      else if (quantity_type == VALUE_DENSITY) {
        string species_name = release_site_item[KEY_MOLECULE].asString();
        if (is_volume_species(data.mcell, species_name)) {
          gen_param_expr(out, NAME_CONCENTRATION, release_site_item[KEY_QUANTITY], prob_not_1);
        }
        else {
          gen_param_expr(out, NAME_DENSITY, release_site_item[KEY_QUANTITY], prob_not_1);
        }
      }
      else {
        ERROR("Quantity type " + quantity_type + " is not supported yet");
      }
    }

    if (prob_not_1) {
      gen_param_expr(out, NAME_RELEASE_PROBABILITY, prob, false);
    }

    out << CTOR_END;
    i++;
  }
}


std::vector<std::string> PythonGenerator::get_species_to_visualize() {
  vector<string> res;

  Value& define_molecules = get_node(mcell, KEY_DEFINE_MOLECULES);
  check_version(KEY_DEFINE_MOLECULES, define_molecules, VER_DM_2014_10_24_1638);

  Value& molecule_list = get_node(define_molecules, KEY_MOLECULE_LIST);
  for (Value::ArrayIndex i = 0; i < molecule_list.size(); i++) {
    Value& molecule_list_item = molecule_list[i];
    check_version(KEY_MOLECULE_LIST, molecule_list_item, VER_DM_2018_10_16_1632);

    if (molecule_list_item[KEY_EXPORT_VIZ].asBool()) {
      if (data.bng_mode) {
        res.push_back(make_species(molecule_list_item[KEY_MOL_NAME].asString()));
      }
      else {
        res.push_back(molecule_list_item[KEY_MOL_NAME].asString());
      }
    }
  }

  return res;
}


void PythonGenerator::generate_viz_outputs(
    std::ostream& out, const bool cellblender_viz,
    std::vector<std::string>& viz_output_names) {

  if (!mcell.isMember(KEY_VIZ_OUTPUT)) {
    return;
  }

  Value& viz_output = mcell[KEY_VIZ_OUTPUT];
  check_version(KEY_VIZ_OUTPUT, viz_output, VER_DM_2014_10_24_1638);

  // species_list
  vector<string> viz_species = get_species_to_visualize();
  if (!viz_output[KEY_EXPORT_ALL].asBool() && viz_species.empty()) {
    return; // nothing to generate
  }

  string name = VIZ_OUTPUT_NAME; // there is only one in datamodel now
  viz_output_names.push_back(name);

  // CHECK_PROPERTY(viz_output[KEY_ALL_ITERATIONS].asBool()); // don't care
  CHECK_PROPERTY(viz_output[KEY_START].asString() == "0");

  data.check_if_already_defined_and_add(name, NAME_CLASS_VIZ_OUTPUT);
  gen_ctor_call(out, name, NAME_CLASS_VIZ_OUTPUT);

  // mode is ascii by default, this information is not in datamodel
  const char* mode = (cellblender_viz) ? NAME_EV_CELLBLENDER : NAME_EV_ASCII;
  gen_param_enum(out, NAME_MODE, NAME_ENUM_VIZ_MODE, mode, true);
  gen_param(out, NAME_OUTPUT_FILES_PREFIX, DEFAULT_VIZ_OUTPUT_FILENAME_PREFIX, true);

  // species_list
  if (!viz_output[KEY_EXPORT_ALL].asBool() && !viz_species.empty()) {
    gen_param_list(out, NAME_SPECIES_LIST, viz_species, true);
  }

  if (viz_output[KEY_STEP].asString() == "0") {
    cout << "Viz output specification has periodicity (" << NAME_EVERY_N_TIMESTEPS << ") value 0. "
        << "This is probably caused by conversion from MDL where the number of simulated iterations "
        << "lower than the actual periodicity. It is not possible to determine the periodicity in this case.\n";
  }
  if (!viz_output[KEY_ALL_ITERATIONS].asBool()) {
    gen_param_expr(out, NAME_EVERY_N_TIMESTEPS, viz_output[KEY_STEP], false);
  }
  else {
    gen_param_expr(out, NAME_EVERY_N_TIMESTEPS, "1", false);
  }

  // ignoring KEY_END
  out << CTOR_END;
}



std::string PythonGenerator::generate_single_count_term(
    ostream& out,
    const std::string& what_to_count,
    const std::string& where_to_count,
    const std::string& orientation,
    const bool molecules_not_species,
    const bool rxn_not_mol
) {
  string name = COUNT_TERM_PREFIX +
      create_count_name(what_to_count, where_to_count, molecules_not_species);

  // generate the count term object definition if we don't already have it
  // TODO: move this into a function
  if (find(data.all_count_term_names.begin(), data.all_count_term_names.end(), name) == data.all_count_term_names.end()) {
    data.all_count_term_names.push_back(name);

    data.check_if_already_defined_and_add(name, NAME_CLASS_COUNT_TERM);

    gen_ctor_call(out, name, NAME_CLASS_COUNT_TERM);

    bool comma_after = (where_to_count != "");

    if (rxn_not_mol) {
      gen_param_expr(out, NAME_REACTION_RULE, what_to_count, comma_after);
    }
    else {
      gen_param_expr(
          out, molecules_not_species ? NAME_MOLECULES_PATTERN : NAME_SPECIES_PATTERN,
          make_species_or_cplx(data, what_to_count, orientation, ""), comma_after);
    }

    if (where_to_count != "") {
      gen_param_expr(out, NAME_REGION, where_to_count, false);
    }

    out << CTOR_END;
  }

  return fix_id(name);
}


// stores multiplier value into the multiplier argument as a string
// if present, the expected form is mult*(<counts>)
string PythonGenerator::generate_count_terms_for_expression(
    ostream& out,
    const string& mdl_string, // may be empty, in that case we use what_to_count and where_to_count
    const std::string& what_to_count,
    const std::string& where_to_count,
    const std::string& orientation,
    const bool rxn_not_mol
) {
  string res_expr;

  if (mdl_string != "") {
    size_t last_end = 0;
    uint num_counts = get_num_counts_in_mdl_string(mdl_string);

    // the first count term item must be positive
    for (uint i = 0; i < num_counts; i++) {
      size_t start = mdl_string.find(COUNT, last_end);

      size_t end = find_end_brace_pos(mdl_string, start + strlen(COUNT));
      if (end == string::npos) {
        end = mdl_string.size();
      }

      bool rxn_not_mol;
      bool molecules_not_species;
      string what_to_count;
      string where_to_count;
      string orientation;

      process_single_count_term(
          data, mdl_string.substr(start, end - start + 1),
          rxn_not_mol, molecules_not_species, what_to_count, where_to_count, orientation
      );

      string name = generate_single_count_term(
          out, what_to_count, where_to_count, orientation, molecules_not_species, rxn_not_mol);

      // for the res_expr, we cut all the COUNT[..] and replace them with
      // the ids of the CountTerm objects
      if (last_end != 0)
        res_expr += " " + mdl_string.substr(last_end + 1, start - (last_end + 1)) + " " + name;
      else {
        res_expr += name;
      }

      last_end = end;
    }
  }
  else {
    bool molecules_not_species = false; // MCell base counting always uses molecules
    string name = COUNT_TERM_PREFIX +
        create_count_name(what_to_count, where_to_count, molecules_not_species);

    res_expr = generate_single_count_term(
        out, what_to_count, where_to_count, orientation, molecules_not_species, rxn_not_mol);
  }

  return res_expr;
}


void PythonGenerator::generate_single_count(
    std::ostream& out,
    const std::string& count_name,
    const std::string& observable_name,
    const std::string& file_name,
    const std::string& count_term_name,
    const std::string& mul_div_str,
    const std::string& rxn_step
) {

  data.check_if_already_defined_and_add(count_name, NAME_CLASS_COUNT);
  gen_ctor_call(out, count_name, NAME_CLASS_COUNT);

  gen_param(out, NAME_NAME, observable_name, true);

  gen_param_expr(out, NAME_EXPRESSION, count_term_name, true);

  gen_param(out, NAME_FILE_NAME,
      DEFAULT_RXN_OUTPUT_FILENAME_PREFIX + file_name, mul_div_str != "" || rxn_step != "");

  if (mul_div_str != "") {
    if (mul_div_str[0] == '*') {
      gen_param_expr(out, NAME_MULTIPLIER, mul_div_str.substr(1), rxn_step != "");
    }
    else if (mul_div_str[0] == '/') {
      gen_param_expr(out, NAME_MULTIPLIER, "(1.0/" + mul_div_str.substr(1) + ")", rxn_step != "");
    }
    else {
      release_assert(false && "Invalid count multiplier or divider.");
    }
  }

  if (rxn_step != "") {
    // rxn_step is specified in seconds, need to convert to seconds
    const string& time_step = mcell[KEY_INITIALIZATION][KEY_TIME_STEP].asString();
    gen_param_expr(out, NAME_EVERY_N_TIMESTEPS, rxn_step + "/" + time_step, false);
  }

  out << CTOR_END;
}


void PythonGenerator::generate_all_bngl_reaction_rules_used_in_observables(std::ostream& out) {

  // bngl_reaction_rules_used_in_observables are initialized in MCell4Generator::generate_reaction_rules
  out << "# ---- declaration of rxn rules defined in BNGL and used in counts ----\n";
  for (string& name: data.bngl_reaction_rules_used_in_observables) {
    out <<
      name << " = " << get_module_name_w_prefix(data.output_files_prefix, SUBSYSTEM) << "." <<
      NAME_FIND_REACTION_RULE << "('" << name << "')\n";
    out << "assert " << name << ", \"Reaction rule '" + name + "' was not found\"\n\n";
  }
  out << "\n";
}


void PythonGenerator::generate_surface_classes_assignments(ostream& out) {
  if (!data.mcell.isMember(KEY_MODIFY_SURFACE_REGIONS)) {
    return;
  }

  Value& modify_surface_regions = data.mcell[KEY_MODIFY_SURFACE_REGIONS];
  check_version(KEY_MODIFY_SURFACE_REGIONS, modify_surface_regions, VER_DM_2014_10_24_1638);
  Value& modify_surface_regions_list = modify_surface_regions[KEY_MODIFY_SURFACE_REGIONS_LIST];
  for (Value::ArrayIndex i = 0; i < modify_surface_regions_list.size(); i++) {
    Value& modify_surface_regions_item = modify_surface_regions_list[i];
    check_versions(
        KEY_MODIFY_SURFACE_REGIONS_LIST, modify_surface_regions_item,
        VER_DM_2018_01_11_1330, VER_DM_2020_07_12_1600
    );

    string object_name = modify_surface_regions_item[KEY_OBJECT_NAME].asString();
    CHECK_PROPERTY(object_name != "");

    string region_selection = modify_surface_regions_item[KEY_REGION_SELECTION].asString();

    string obj_or_region_name;
    if (region_selection == VALUE_ALL) {
      obj_or_region_name = object_name;
    }
    else if (region_selection == VALUE_SEL) {
      obj_or_region_name = object_name + "_" + modify_surface_regions_item[KEY_REGION_NAME].asString();
    }
    else {
      ERROR("Unexpected value " + region_selection + " of " + KEY_REGION_SELECTION + ".");
    }

    // handle initial_region_molecules_list if present
    if (modify_surface_regions_item.isMember(KEY_INITIAL_REGION_MOLECULES_LIST)) {
      string initial_region_molecules_name = obj_or_region_name + "_" + NAME_INITIAL_SURFACE_RELEASES;
      out << initial_region_molecules_name << " = [\n";

      Value& initial_region_molecules_list = modify_surface_regions_item[KEY_INITIAL_REGION_MOLECULES_LIST];
      for (Value::ArrayIndex rel_i = 0; rel_i < initial_region_molecules_list.size(); rel_i++) {
        Value& item = initial_region_molecules_list[rel_i];
        out << "    ";
        gen_ctor_call(out, "", NAME_CLASS_INITIAL_SURFACE_RELEASE, true);
        out << "    ";
        string species_str = item[KEY_MOLECULE].asString();
        string orient = convert_orientation(item[KEY_ORIENT].asString(), true);
        gen_param_expr(out, NAME_COMPLEX, make_species_or_cplx(data, species_str, orient), true);
        out << "    ";
        if (item.isMember(KEY_MOLECULE_NUMBER)) {
          gen_param_expr(out, NAME_NUMBER_TO_RELEASE, item[KEY_MOLECULE_NUMBER], false);
        }
        else if (item.isMember(KEY_MOLECULE_DENSITY)) {
          gen_param_expr(out, NAME_DENSITY, item[KEY_MOLECULE_DENSITY], false);
        }
        else {
          ERROR(
              S("Missing ") + KEY_MOLECULE_NUMBER + " or " + KEY_MOLECULE_DENSITY +
              " in " + KEY_INITIAL_REGION_MOLECULES_LIST + "."
          );
        }
        out << "    )";
        gen_comma(out, rel_i, initial_region_molecules_list);
        out << "\n";
      }
      out << "]\n\n";

      out << obj_or_region_name << "." << NAME_INITIAL_SURFACE_RELEASES << " = " << initial_region_molecules_name << "\n";
    }

    // and also surface class name
    if (modify_surface_regions_item.isMember(KEY_SURF_CLASS_NAME)) {

      gen_description(out, modify_surface_regions_item);

      string surf_class_name = modify_surface_regions_item[KEY_SURF_CLASS_NAME].asString();
      out << obj_or_region_name << "." << NAME_SURFACE_CLASS << " = " << surf_class_name << "\n";
    }
  }
}


void PythonGenerator::generate_compartment_assignments(std::ostream& out) {
  Value& model_objects = get_node(data.mcell, KEY_MODEL_OBJECTS);
  Value& model_object_list = get_node(model_objects, KEY_MODEL_OBJECT_LIST);
  if (model_object_list.empty()) {
    return;
  }

  for (Value::ArrayIndex i = 0; i < model_object_list.size(); i++) {
    Value& model_object = model_object_list[i];

    // older data model files don't have to have this attribute
    // therefore we are also checking for used compartments in
    // data.is_used_compartment
    bool is_bngl_compartment =
        model_object.isMember(KEY_IS_BNGL_COMPARTMENT) && model_object[KEY_IS_BNGL_COMPARTMENT].asBool();

    if (is_bngl_compartment || data.is_used_compartment(model_object)) {
      // name is the name of the object
      const string& name = model_object[KEY_NAME].asString();
      const string& membrane_name = model_object[KEY_MEMBRANE_NAME].asString();

      gen_assign(out, name, NAME_IS_BNGL_COMPARTMENT, true);
      if (membrane_name != "") {
        gen_assign_str(out, name, NAME_SURFACE_COMPARTMENT_NAME, membrane_name);
        data.surface_to_volume_compartments_map[membrane_name] = name;
      }
      out << "\n";
    }
  }
}

} /* namespace MCell */
