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

#include <algorithm>

#include "generator_utils.h"
#include "mcell4_generator.h"

using namespace std;
using namespace MCell::API;

namespace MCell {

using Json::Value;


string MCell4Generator::get_module_name(const string file_suffix) {
  return get_module_name_w_prefix(output_files_prefix, file_suffix);
}


string MCell4Generator::make_import(const string file_suffix) {
  return "from " + get_module_name(file_suffix) + " import *\n";
}


void MCell4Generator::open_and_check_file(
    const string file_suffix, ofstream& out,
    const bool for_append, const bool bngl) {

  open_and_check_file_w_prefix(output_files_prefix, file_suffix, out, for_append, bngl);
}


void MCell4Generator::reset() {
  unnamed_rxn_counter = 0;
  geometry_generated = false;
  observables_generated = false;
  all_species_and_mol_type_names.clear();
  all_reaction_rules_names.clear();
  bngl_reaction_rules_used_in_observables.clear();
  all_count_term_names.clear();
  bng_gen = nullptr;
  python_gen = nullptr;
}


// the aim is to generate as much of output as possible,
// therefore we are using exceptions
bool MCell4Generator::generate(
    const string& input_file,
    const string& output_files_prefix_,
    const bool bng_mode_,
    const bool debug_mode_,
    const bool cellblender_viz
) {
  reset();

  bool failed = false;
  output_files_prefix = output_files_prefix_;
  bng_mode = bng_mode_;
  debug_mode = debug_mode_;

  // load json file
  ifstream file;
  file.open(input_file);
  if (!file.is_open()) {
    cerr << "Could not open file '" << input_file << "' for reading.\n";
    return false;
  }
  Value root;
  file >> root;
  file.close();

  mcell = get_node(KEY_ROOT, root, KEY_MCELL);

  // create generators
  if (bng_mode) {
    open_and_check_file(MODEL, bng_out, false, true);
    bng_gen = new BNGLGenerator(
        get_filename(output_files_prefix, MODEL, BNGL_EXT), bng_out,
        mcell, output_files_prefix, unnamed_rxn_counter);
  }
  python_gen = new PythonGenerator(mcell, output_files_prefix, bng_mode, unnamed_rxn_counter);

  CHECK(check_scripting(), failed);

  CHECK(generate_parameters(), failed);
  CHECK(generate_subsystem(), failed);
  std::vector<std::string> geometry_names;
  CHECK(geometry_names = generate_geometry(), failed);
  CHECK(generate_instantiation(geometry_names), failed);
  CHECK(generate_observables(cellblender_viz), failed);
  CHECK(generate_model(failed), failed);

  // delete generators
  if (bng_mode) {
    delete bng_gen;
    bng_gen = nullptr;
    bng_out.close();
  }
  delete python_gen;
  python_gen = nullptr;

  return !failed;
}

void MCell4Generator::check_scripting() {
  if (mcell.isMember(KEY_SCRIPTING)) {
    Value& scripting = get_node(mcell, KEY_SCRIPTING);
    if (scripting.isMember(KEY_SCRIPTING_LIST)) {
      Value& scripting_list = get_node(scripting, KEY_SCRIPTING_LIST);
      if (!scripting_list.empty()) {
        ERROR("Data model contains scripting. To convert this model, generate first MDL, "
            "then  convert MDL to data model and then use this converter. "
            "Conversion will continue, but it will be wrong."
        );
      }
    }
  }
}


void MCell4Generator::generate_parameters() {
  ofstream out;
  open_and_check_file(PARAMETERS, out);

  out << make_section_comment("model parameters");

  if (mcell.isMember(KEY_PARAMETER_SYSTEM)) {
    if (bng_mode) {
      bng_gen->generate_parameters(out);
    }
    else {
      python_gen->generate_parameters(out);
    }
  }

  out << make_section_comment("simulation setup");

  out << PARAM_ITERATIONS << " = " << mcell[KEY_INITIALIZATION][KEY_ITERATIONS].asString() << "\n";
  out << PARAM_TIME_STEP << " = " << mcell[KEY_INITIALIZATION][KEY_TIME_STEP].asString() << "\n";
  out << PARAM_DUMP << " = " << (debug_mode ? "True" : "False") << "\n";
  out << PARAM_EXPORT_DATA_MODEL << " = " << "True\n";
  out << "\n";

  out <<
      "# do not use the variable " PARAM_SEED " directly,\n"
      "# Python on import creates copies that do not reflect the current value\n"
      PARAM_SEED << " = 1\n"
      "\n"
      "def " FUNCTION_UPDATE_SEED "(new_value):\n"
      "    global " PARAM_SEED "\n"
      "    " PARAM_SEED " = new_value\n"
      "\n"
      "def " FUNCTION_GET_SEED "():\n"
      "    return " PARAM_SEED "\n"
      "\n";

  out.close();
}


void MCell4Generator::generate_species_and_mol_types(
    ofstream& out, vector<SpeciesOrMolType>& species_and_mt_info) {

  // skip if there are no species/mol types
  if (!mcell.isMember(KEY_DEFINE_MOLECULES)) {
    return;
  }

  if (bng_mode) {
    // molecule types are optional but they allow for better BNGL semantic checks

    // - we also need to generate code that sets diffusion constant,
    //   custom time/space step and other parameters to molecule types,
    // - it must be executed after bngl file is loaded, so we should not be putting it into
    //   subsystem, the new file can be called bngl_molecule_types_info.py
    // - all this extra information should be put preferably into the BNGL file
    //   but for now we would like it to be BNGL compatible and parsing comments or MCELL_ parameters is not really nice
    ofstream mt_info_out;
    open_and_check_file(BNGL_MOLECULE_TYPES_INFO, mt_info_out);
    mt_info_out << MCELL_IMPORT;
    mt_info_out << make_import(PARAMETERS);
    mt_info_out << "\n";

    bng_gen->generate_mol_types(mt_info_out);

    mt_info_out.close();
  }
  else {
    python_gen->generate_species_and_mol_types(out, species_and_mt_info);
  }
}


static bool rxn_uses_mcell_orientation(Value& reaction_list_item) {
  string substances =
      reaction_list_item[KEY_REACTANTS].asString() + " " +
      reaction_list_item[KEY_PRODUCTS].asString();

  // looking for ',; outside of parentheses
  bool in_paren = false;
  size_t i = 0;
  while (i < substances.size()) {
    char c = substances[i];
    if (c == '(') {
      CHECK_PROPERTY(!in_paren && "Malformed reaction definition - embedded parentheses");
      in_paren = true;
    }
    else if (c == ')') {
      CHECK_PROPERTY(in_paren && "Malformed reaction definition - unexpected closing parenthesis");
      in_paren = false;
    }
    else if (!in_paren && (c == ',' || c == '\'' || c == ';')) {
      return true;
    }

    i++;
  }

  return false;
}


static bool rxn_has_variable_rate(Value& reaction_list_item) {
  return reaction_list_item[KEY_VARIABLE_RATE_SWITCH].asBool();
}


// returns true if the rxn name might be referenced by counts
static bool is_rxn_used_in_observables(Value& mcell, const string& rxn_name) {
  Value& reaction_data_output = get_node(mcell, KEY_REACTION_DATA_OUTPUT);
  Value& reaction_output_list = get_node(reaction_data_output, KEY_REACTION_OUTPUT_LIST);
  for (Value::ArrayIndex i = 0; i < reaction_output_list.size(); i++) {
    Value& reaction_output_item = reaction_output_list[i];
    string count_mdl_string = reaction_output_item[KEY_MDL_STRING].asString();
    string count_rxn_name = reaction_output_item[KEY_REACTION_NAME].asString();
    string string_to_check = count_mdl_string + " " +count_rxn_name;

    if (string_to_check.find(rxn_name) != string::npos) {
      return true;
    }
  }
  return false;
}


vector<IdLoc> MCell4Generator::generate_reaction_rules(ofstream& out) {
  vector<IdLoc> rxn_names_w_loc;

  if (!mcell.isMember(KEY_DEFINE_REACTIONS)) {
    return rxn_names_w_loc;
  }

  if (bng_mode) {
    // put into BNG all rxn rules that,
    // 1) don't use MCell orientation or
    // 2) variable rxn rates
    // rest goes into python

    Value& define_reactions = get_node(mcell, KEY_DEFINE_REACTIONS);
    check_version(KEY_DEFINE_REACTIONS, define_reactions, VER_DM_2014_10_24_1638);

    bng_gen->open_reaction_rules_section();

    Value& reaction_list = get_node(define_reactions, KEY_REACTION_LIST);
    for (Value::ArrayIndex i = 0; i < reaction_list.size(); i++) {
      Value& reaction_list_item = reaction_list[i];
      check_version(KEY_MOLECULE_LIST, reaction_list_item, VER_DM_2018_01_11_1330);

      if (!rxn_uses_mcell_orientation(reaction_list_item) && !rxn_has_variable_rate(reaction_list_item)) {

        bool used_in_observables = is_rxn_used_in_observables(mcell, reaction_list_item[KEY_RXN_NAME].asString());
        string name = bng_gen->generate_single_reaction_rule(reaction_list_item, used_in_observables);
        rxn_names_w_loc.push_back(IdLoc(name, false));
        bngl_reaction_rules_used_in_observables.push_back(name);
      }
      else {
        bng_gen->add_comment(
            S(IND) + "reaction " + reaction_list_item[KEY_MOL_NAME].asString() +
            " was generated as Python code because it contains features not supported by BNGL");
        string name = python_gen->generate_single_reaction_rule(out, reaction_list_item);
        rxn_names_w_loc.push_back(IdLoc(name, true));
      }
    }

    bng_gen->close_reaction_rules_section();
  }
  else {
    python_gen->generate_reaction_rules(out, rxn_names_w_loc);
  }

  return rxn_names_w_loc;
}


void MCell4Generator::generate_subsystem() {
  ofstream out;
  open_and_check_file(SUBSYSTEM, out);

  out << MCELL_IMPORT;
  out << make_import(PARAMETERS);
  out << "\n";
  out << make_section_comment(SUBSYSTEM);

  generate_species_and_mol_types(out, all_species_and_mol_type_names);

  vector<string> surface_class_names;
  python_gen->generate_surface_classes(out, surface_class_names);

  all_reaction_rules_names = generate_reaction_rules(out);

  gen_ctor_call(out, SUBSYSTEM, NAME_CLASS_SUBSYSTEM, false);
  for (SpeciesOrMolType& info: all_species_and_mol_type_names) {
    if (info.is_species) {
      gen_method_call(out, SUBSYSTEM, NAME_ADD_SPECIES, info.name);
    }
    else {
      gen_method_call(out, SUBSYSTEM, NAME_ADD_ELEMENTARY_MOLECULE_TYPE, info.name);
    }
  }
  for (string& sc: surface_class_names) {
    gen_method_call(out, SUBSYSTEM, NAME_ADD_SURFACE_CLASS, sc);
  }
  for (IdLoc& r_loc: all_reaction_rules_names) {
    if (r_loc.in_python) {
      gen_method_call(out, SUBSYSTEM, NAME_ADD_REACTION_RULE, r_loc.name);
    }
  }

  out.close();
}


vector<string> MCell4Generator::generate_geometry() {

  vector<string> geometry_objects;

  // TODO: check versions

  if (!mcell.isMember(KEY_GEOMETRICAL_OBJECTS)) {
    return geometry_objects;
  }
  Value& geometrical_objects = get_node(mcell, KEY_GEOMETRICAL_OBJECTS);
  if (!geometrical_objects.isMember(KEY_OBJECT_LIST)) {
    return geometry_objects;
  }

  ofstream out;
  open_and_check_file(GEOMETRY, out);

  out << MCELL_IMPORT;

  python_gen->generate_geometry(out, geometry_objects);

  // NOTE: we can generate BNGL compartments from geometry here

  out.close();
  geometry_generated = true;

  return geometry_objects;
}



void MCell4Generator::generate_surface_classes_assignment(ofstream& out) {
  if (!mcell.isMember(KEY_MODIFY_SURFACE_REGIONS)) {
    return;
  }

  Value& modify_surface_regions = mcell[KEY_MODIFY_SURFACE_REGIONS];
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
        gen_param_id(out, NAME_SPECIES, item[KEY_MOLECULE], true);
        out << "    ";
        gen_param_enum(out, NAME_ORIENTATION, NAME_ENUM_ORIENTATION, convert_orientation(item[KEY_ORIENT].asString()), true);
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
      string surf_class_name = modify_surface_regions_item[KEY_SURF_CLASS_NAME].asString();
      out << obj_or_region_name << "." << NAME_SURFACE_CLASS << " = " << surf_class_name << "\n";
    }

  }
  out << "\n";
}


void MCell4Generator::generate_instantiation(const vector<string>& geometry_objects) {

  ofstream out;
  open_and_check_file(INSTANTIATION, out);

  out << MCELL_IMPORT;
  out << make_import(PARAMETERS);
  out << make_import(SUBSYSTEM);
  if (geometry_generated) {
    out << make_import(GEOMETRY);
  }
  out << "\n";
  out << make_section_comment(INSTANTIATION);
  out << make_section_comment("release sites");

  // BNGL instantiation through seed species are not supported yet, everything is defined through Python

  vector<string> release_sites;
  python_gen->generate_release_sites(out, release_sites);

  out << make_section_comment("surface classes assignment");
  generate_surface_classes_assignment(out);

  out << make_section_comment("instantiation data");

  gen_ctor_call(out, INSTANTIATION, NAME_CLASS_INSTANTIATION_DATA, false);
  for (const string& s: geometry_objects) {
    gen_method_call(out, INSTANTIATION, NAME_ADD_GEOMETRY_OBJECT, s);
  }
  for (const string& r: release_sites) {
    gen_method_call(out, INSTANTIATION, NAME_ADD_RELEASE_SITE, r);
  }

  out.close();
}


static uint get_num_counts_in_mdl_string(const string& mdl_string) {
  uint res = 0;
  size_t pos = 0;
  while ((pos = mdl_string.find(COUNT, pos)) != string::npos) {
    res++;
    pos += strlen(COUNT);
  }
  return res;
}


static string create_count_name(string what_to_count, string where_to_count) {

  // first remove all cterm_refixes
  regex pattern_cterm(COUNT_TERM_PREFIX);
  what_to_count = regex_replace(what_to_count, pattern_cterm, "");


  string res = COUNT_PREFIX;
  for (char c: what_to_count) {
    if (c == '+') {
      res += "_plus_";
    }
    else if (c == '-') {
      res += "_minus_";
    }
    else if (c == '(') {
      res += "_pstart_";
    }
    else if (c == ')') {
      res += "_pend_";
    }
    else if (c == '.') {
      res += "_";
    }
    else if (c == '_') {
      res += "_";
    }
    else if (isalnum(c)) {
      res += c;
    }
    // ignoring the rest of the characters
  }

  if (where_to_count != WORLD && where_to_count != "") {
    res += "_" + where_to_count;
  }

  return res;
}


void MCell4Generator::process_single_count_term(
    const string& mdl_string,
    bool& rxn_not_mol, string& what_to_count, string& where_to_count, string& orientation) {

  // mdl_string is always in the form COUNT[what,where]
  size_t start_brace = mdl_string.find('[');
  size_t comma = mdl_string.find(',');
  size_t comma2 = mdl_string.find(',', comma + 1);
  if (comma2 != string::npos) {
    // the first comma is orientation
    comma = comma2;
  }
  size_t end_brace = mdl_string.find(']');

  if (mdl_string.find(COUNT) == string::npos) {
    ERROR("String 'COUNT' was not found in mdl_string '" + mdl_string + "'.");
  }
  if (start_brace == string::npos || comma == string::npos || end_brace == string::npos ||
      start_brace > comma || start_brace > end_brace || comma > end_brace
  ) {
    ERROR("Malformed mdl_string '" + mdl_string + "'.");
  }

  what_to_count = mdl_string.substr(start_brace + 1, comma - start_brace - 1);
  what_to_count = trim(what_to_count);

  // process orientation
  orientation = "";
  assert(what_to_count != "");
  char last_c = what_to_count.back();
  if (last_c == '\'' || last_c == ',' || last_c == ';') {
    string s;
    s = last_c;
    orientation = convert_orientation(s);

    what_to_count = what_to_count.substr(0, what_to_count.size() - 1);
  }

  if (find(all_species_and_mol_type_names.begin(), all_species_and_mol_type_names.end(), SpeciesOrMolType(what_to_count))
      != all_species_and_mol_type_names.end()) {
    rxn_not_mol = false;
  }
  else if (find(all_reaction_rules_names.begin(), all_reaction_rules_names.end(), IdLoc(what_to_count)) != all_reaction_rules_names.end()) {
    rxn_not_mol = true;
  }
  else {
    ERROR("Identifier '" + what_to_count + "' is neither species nor reaction, from mdl_string '" + mdl_string + "'.");
  }

  where_to_count = mdl_string.substr(comma + 1, end_brace - comma - 1);
  size_t dot_pos = where_to_count.find('.');
  if (dot_pos !=- string::npos) {
    where_to_count = where_to_count.substr(dot_pos + 1);
  }
  where_to_count = trim(where_to_count);
  if (where_to_count == WORLD) {
    where_to_count = "";
  }
  // where_to_count can now look like this: "Cube[ALL"
  size_t brace = where_to_count.find('[');
  if (brace != string::npos) {
    if (where_to_count.substr(brace + 1) == REGION_ALL_NAME) {
      where_to_count = where_to_count.substr(0, brace);
    }
    else {
      where_to_count[brace] = '_';
    }
  }
}


static string get_count_multiplier(const string& mdl_string) {
  // check if there is a multiplier
  // handling only code generated by mcell4 for now
  string multiplier_str = "";
  size_t mult_pos = mdl_string.find("*");
  if (mult_pos != string::npos) {
    string substr = mdl_string.substr(0, mult_pos);
    try {
      stod(substr);
      multiplier_str = substr;
    }
    catch (const std::exception& e) {
      ERROR("Could not convert multiplier from " + mdl_string + ".");
    }
  }
  return multiplier_str;
}

// stores multiplier value into the multiplier argument as a string
// if present, the expected form is mult*(<counts>)
string MCell4Generator::generate_count_terms_for_expression(
    ofstream& out, const string& mdl_string) {
  string res_expr;

  size_t last_end = 0;
  uint num_counts = get_num_counts_in_mdl_string(mdl_string);

  // the first count term item must be positive
  for (uint i = 0; i < num_counts; i++) {
    size_t start = mdl_string.find(COUNT, last_end);

    size_t end = mdl_string.find(']', start);
    if (end == string::npos) {
      end = mdl_string.size();
    }

    bool rxn_not_mol;
    string what_to_count;
    string where_to_count;
    string orientation;

    process_single_count_term(
        mdl_string.substr(start, end - start + 1),
        rxn_not_mol, what_to_count, where_to_count, orientation
    );

    string name = COUNT_TERM_PREFIX + what_to_count + ((where_to_count != "") ? ("_" + where_to_count) : "");

    // generate the count term object definition if we don't already have it
    if (find(all_count_term_names.begin(), all_count_term_names.end(), name) == all_count_term_names.end()) {
      all_count_term_names.push_back(name);
      gen_ctor_call(out, name, NAME_CLASS_COUNT_TERM);

      if (rxn_not_mol) {
        gen_param_id(out, NAME_REACTION_RULE, what_to_count, where_to_count != "");
      }
      else {
        gen_param_id(out, NAME_SPECIES, what_to_count, orientation == "" && where_to_count != "");

        if (orientation != "") {
          gen_param_enum(out, NAME_ORIENTATION, NAME_ENUM_ORIENTATION, orientation, where_to_count != "");
        }
      }

      if (where_to_count != "") {
        gen_param_id(out, NAME_REGION, where_to_count, false);
      }

      out << CTOR_END;
    }

    // for the res_expr, we cut all the COUNT[..] and replace them with
    // the ids of the CountTerm objects
    if (last_end != 0)
      res_expr += " " + mdl_string.substr(last_end + 1, start - (last_end + 1)) + " " + name;
    else {
      res_expr += name;
    }

    last_end = end;
  }

  return res_expr;
}


vector<string> MCell4Generator::generate_counts(ofstream& out) {
  vector<string> counts;

  if (!mcell.isMember(KEY_REACTION_DATA_OUTPUT)) {
    return counts;
  }

  Value& reaction_data_output = get_node(mcell, KEY_REACTION_DATA_OUTPUT);
  check_version(KEY_DEFINE_MOLECULES, reaction_data_output, VER_DM_2016_03_15_1800);

  string rxn_step = reaction_data_output[KEY_RXN_STEP].asString();

  Value& reaction_output_list = get_node(reaction_data_output, KEY_REACTION_OUTPUT_LIST);
  for (Value::ArrayIndex i = 0; i < reaction_output_list.size(); i++) {
    Value& reaction_output_item = reaction_output_list[i];

    string mdl_file_prefix = reaction_output_item[KEY_MDL_FILE_PREFIX].asString();

    bool rxn_not_mol;
    bool single_term;
    string what_to_count;
    string where_to_count; // empty for WORLD
    string orientation;

    string count_location = reaction_output_item[KEY_COUNT_LOCATION].asString();

    string rxn_or_mol = reaction_output_item[KEY_RXN_OR_MOL].asString();
    string multiplier_str = "";
    if (rxn_or_mol == VALUE_MDLSTRING) {
      // first check whether we need to generate count_terms
      string mdl_string = reaction_output_item[KEY_MDL_STRING].asString();
      uint num_counts = get_num_counts_in_mdl_string(mdl_string);
      if (num_counts == 0) {
        ERROR("There is no 'COUNT' in mdl_string for output with filename " +  mdl_file_prefix + ".");
      }
      else if (num_counts == 1) {
        single_term = true;
        process_single_count_term(mdl_string, rxn_not_mol, what_to_count, where_to_count, orientation);
      }
      else {
        single_term = false;
        what_to_count = generate_count_terms_for_expression(out, mdl_string);
        where_to_count = "";
      }

      multiplier_str = get_count_multiplier(mdl_string);
    }
    else if (rxn_or_mol == VALUE_REACTION) {
      single_term = true;
      rxn_not_mol = true;
      what_to_count = reaction_output_item[KEY_REACTION_NAME].asString();

      if (count_location == VALUE_COUNT_LOCATION_OBJECT) {
        where_to_count = reaction_output_item[KEY_OBJECT_NAME].asString();
      }
      else if (count_location == VALUE_COUNT_LOCATION_REGION) {
        where_to_count = reaction_output_item[KEY_REGION_NAME].asString();
      }
      else {
        assert(count_location == VALUE_COUNT_LOCATION_WORLD);
      }
    }
    else if (rxn_or_mol == VALUE_MOLECULE) {
      single_term = true;
      rxn_not_mol = false;
      what_to_count = reaction_output_item[KEY_MOLECULE_NAME].asString();

      if (count_location == VALUE_COUNT_LOCATION_OBJECT) {
        where_to_count = reaction_output_item[KEY_OBJECT_NAME].asString();
      }
      else if (count_location == VALUE_COUNT_LOCATION_REGION) {
        where_to_count = reaction_output_item[KEY_REGION_NAME].asString();
      }
      else {
        assert(count_location == VALUE_COUNT_LOCATION_WORLD);
      }
    }
    else {
      ERROR("Invalid rxn_or_mol '" + rxn_or_mol + "' in reaction_output_list for output with filename " +
          mdl_file_prefix + ".");
    }

    string name = create_count_name(what_to_count, where_to_count);
    counts.push_back(name);
    gen_ctor_call(out, name, NAME_CLASS_COUNT);

    if (single_term) {
      if (rxn_not_mol) {
        gen_param_id(out, NAME_REACTION_RULE, what_to_count, true);
      }
      else {
        gen_param_id(out, NAME_SPECIES, what_to_count, true);

        if (orientation != "") {
          gen_param_enum(out, NAME_ORIENTATION, NAME_ENUM_ORIENTATION, orientation, where_to_count != "");
        }
      }
    }
    else {
      gen_param_id(out, NAME_COUNT_EXPRESSION, what_to_count, true);
    }

    if (where_to_count != "") {
      gen_param_id(out, NAME_REGION, where_to_count, true);
    }

    if (mdl_file_prefix == "") {
      string where = where_to_count;
      if (where == "") {
        where = WORLD_FIRST_UPPER;
      }
      mdl_file_prefix = what_to_count + "." + where;

      // TODO: this might need further checks
      if (mdl_file_prefix.find_first_of(" ,+*/\\") != string::npos) {
        cout << "Warning: count file prefix '" + mdl_file_prefix + "' is probably invalid.\n";
      }
    }

    gen_param(out, NAME_FILE_NAME,
        DEFAULT_RXN_OUTPUT_FILENAME_PREFIX + mdl_file_prefix + ".dat", multiplier_str != "" || rxn_step != "");

    if (multiplier_str != "") {
      gen_param_expr(out, NAME_MULTIPLIER, multiplier_str, rxn_step != "");
    }

    if (rxn_step != "") {
      gen_param_expr(out, NAME_EVERY_N_TIMESTEPS, rxn_step, false);
    }

    out << CTOR_END;
  }

  return counts;
}


void MCell4Generator::generate_observables(const bool cellblender_viz) {

  ofstream out;
  open_and_check_file(OBSERVABLES, out);

  out << MCELL_IMPORT;
  out << make_import(PARAMETERS);
  out << make_import(SUBSYSTEM);
  if (geometry_generated) {
    out << make_import(GEOMETRY);
  }
  out << "\n";
  out << make_section_comment(OBSERVABLES);

  vector<string> viz_outputs;
  python_gen->generate_viz_outputs(out, cellblender_viz, viz_outputs);

  vector<string> counts = generate_counts(out);

  gen_ctor_call(out, OBSERVABLES, NAME_CLASS_OBSERVABLES, false);
  for (const string& s: viz_outputs) {
    gen_method_call(out, OBSERVABLES, NAME_ADD_VIZ_OUTPUT, s);
  }
  for (const string& r: counts) {
    gen_method_call(out, OBSERVABLES, NAME_ADD_COUNT, r);
  }

  observables_generated = true;

  out.close();
}


void MCell4Generator::generate_config(ofstream& out) {
  out << make_section_comment("configuration");

  // using values from generated parameters.py
  gen_assign(out, MODEL, NAME_CONFIG, NAME_TIME_STEP, PARAM_TIME_STEP);
  gen_assign(out, MODEL, NAME_CONFIG, NAME_SEED, S(FUNCTION_GET_SEED) + "()");
  gen_assign(out, MODEL, NAME_CONFIG, NAME_TOTAL_ITERATIONS_HINT, PARAM_ITERATIONS);
  out << "\n";

  if (!mcell.isMember(KEY_INITIALIZATION)) {
    ERROR(S("Data model does not contain key ") + KEY_INITIALIZATION + ".");
  }

  Value& initialization = mcell[KEY_INITIALIZATION];
  Value& partitions = initialization[KEY_PARTITIONS];

  // choose the largest value for partition size and the smallest step
  float_t x_start = stod(partitions[KEY_X_START].asString());
  float_t x_end = stod(partitions[KEY_X_END].asString());
  float_t x_step = stod(partitions[KEY_X_STEP].asString());
  float_t y_start = stod(partitions[KEY_Y_START].asString());
  float_t y_end = stod(partitions[KEY_Y_END].asString());
  float_t y_step = stod(partitions[KEY_Y_STEP].asString());
  float_t z_start = stod(partitions[KEY_Z_START].asString());
  float_t z_end = stod(partitions[KEY_Z_END].asString());
  float_t z_step = stod(partitions[KEY_Z_STEP].asString());

  float_t partition_dimension;
  if (!cmp_eq(x_start, -x_end) || !cmp_eq(y_start, -y_end) || !cmp_eq(y_start, -y_end) ||
      !cmp_eq(x_start, y_start) || !cmp_eq(y_start, z_start) ||
      !cmp_eq(x_end, y_end) || !cmp_eq(y_end, z_end)
  ) {
    gen_assign_vec3(out, MODEL, NAME_CONFIG, NAME_INITIAL_PARTITION_ORIGIN, x_start, y_start, z_start);
    // select the largest difference as partition dimension
    Vec3 dims = Vec3(x_end - x_start, y_end - y_start, z_end - z_start);
    partition_dimension = max3(dims);
  }
  else {
    partition_dimension = x_end * 2;
  }


  float_t partition_step;
  if (!cmp_eq(x_step, y_step) || !cmp_eq(y_step, z_step)) {
    cout << "Message: Partition's step sizes are different, changing the step to be the smallest of them.";

    partition_step = min3d(x_step, y_step, z_step);
  }
  else {
    partition_step = x_step;
  }

  gen_assign(out, MODEL, NAME_CONFIG, NAME_PARTITION_DIMENSION, partition_dimension);
  gen_assign(out, MODEL, NAME_CONFIG, NAME_SUBPARTITION_DIMENSION, partition_step);

  out << "\n";
  out << make_section_comment("default configuration overrides");

  // FIXME: check the value against the default in the yaml API file
  bool center_molecules_on_grid = initialization[KEY_CENTER_MOLECULES_ON_GRID].asBool();
  if (center_molecules_on_grid) {
    gen_assign(out, MODEL, NAME_CONFIG, NAME_CENTER_MOLECULES_ON_GRID, center_molecules_on_grid);
  }

  string radius = initialization[KEY_INTERACTION_RADIUS].asString();
  if (radius != "") {
    gen_assign(out, MODEL, NAME_CONFIG, NAME_CENTER_MOLECULES_ON_GRID, stod(radius));
  }

  string surf_grid_density = initialization[KEY_SURFACE_GRID_DENSITY].asString();
  if (surf_grid_density != "" && surf_grid_density != "10000") {
    gen_assign(out, MODEL, NAME_CONFIG, NAME_CENTER_MOLECULES_ON_GRID, stod(surf_grid_density));
  }

  string vacancy_search_distance = initialization[KEY_VACANCY_SEARCH_DISTANCE].asString();
  if (vacancy_search_distance != "" && vacancy_search_distance != "10") {
    gen_assign(out, MODEL, NAME_CONFIG, NAME_VACANCY_SEARCH_DISTANCE, stod(vacancy_search_distance));
  }
}


void MCell4Generator::generate_model(const bool print_failed_marker) {
  ofstream out;
  open_and_check_file(MODEL, out);

  out << INTERPRETER;

  if (print_failed_marker) {
    out <<
        "ERROR: Conversion from data model failed, these generated sources are result of a "
        "best-effort conversion and will contain errors that might need to be fixed manually.\n\n";
  }

  out << BASE_MODEL_IMPORTS;
  out << MCELL_DIR_SETUP;

  out << MCELL_IMPORT;

  out << make_import(PARAMETERS);

  out <<
      "\n"
      "if len(sys.argv) == 3 and sys.argv[1] == '-seed\':\n"
      "    # overwrite value of seed defined in module " + get_module_name(PARAMETERS) + "\n"
      "    " FUNCTION_UPDATE_SEED "(int(sys.argv[2]))\n\n";

  out << IMPORT << " " << get_module_name(SUBSYSTEM) << "\n";
  out << IMPORT << " " << get_module_name(INSTANTIATION) << "\n";
  if (observables_generated) {
    out << IMPORT << " " << get_module_name(OBSERVABLES) << "\n";
  }
  out << "\n";

  gen_ctor_call(out, MODEL, NAME_CLASS_MODEL, false);
  out << "\n";

  generate_config(out);
  out << "\n";

  out << make_section_comment("add components");
  gen_method_call(out, MODEL, NAME_ADD_SUBSYSTEM, get_module_name(SUBSYSTEM) + "." + SUBSYSTEM);
  gen_method_call(out, MODEL, NAME_ADD_INSTANTIATION_DATA, get_module_name(INSTANTIATION) + "." + INSTANTIATION);
  if (observables_generated) {
    gen_method_call(out, MODEL, NAME_ADD_OBSERVABLES, get_module_name(OBSERVABLES) + "." + OBSERVABLES);
  }
  out << "\n";

  out << make_section_comment("initialization and execution");
  gen_method_call(out, MODEL, NAME_INITIALIZE);
  out << "\n";

  out << "if " << PARAM_DUMP << ":\n";
  out << IND;
  gen_method_call(out, MODEL, NAME_DUMP_INTERNAL_STATE);
  out << "\n";

  // method export_data_model uses target directory from viz_outputs
  out << "if " << PARAM_EXPORT_DATA_MODEL << " and " << MODEL << "." << NAME_VIZ_OUTPUTS << ":\n";
  out << IND;
  gen_method_call(out, MODEL, NAME_EXPORT_DATA_MODEL);
  out << "\n";

  gen_method_call(out, MODEL, NAME_RUN_ITERATIONS, PARAM_ITERATIONS);
  gen_method_call(out, MODEL, NAME_END_SIMULATION);
}


} // namespace MCell
