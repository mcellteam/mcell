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
#include "bng/bng_defines.h"
#include "src/util.h"

using namespace std;
using namespace MCell::API;

namespace MCell {

using Json::Value;


string MCell4Generator::get_module_name(const string file_suffix) {
  return get_module_name_w_prefix(data.output_files_prefix, file_suffix);
}


string MCell4Generator::make_import(const string file_suffix) {
  return "from " + get_module_name(file_suffix) + " import *\n";
}


void MCell4Generator::open_and_check_file(
    const string file_suffix, ofstream& out,
    const bool for_append, const bool bngl) {

  open_and_check_file_w_prefix(data.output_files_prefix, file_suffix, out, for_append, bngl);
}


void MCell4Generator::reset() {
  data.reset();
  geometry_generated = false;
  observables_generated = false;
  bng_gen = nullptr;
  python_gen = nullptr;
}


// the aim is to generate as much of output as possible,
// therefore we are using exceptions
bool MCell4Generator::generate(const SharedGenData& opts) {
  reset();

  bool failed = false;
  data = opts; // copy options

  // load json file
  ifstream file;
  file.open(opts.input_file);
  if (!file.is_open()) {
    cerr << "Could not open file '" << opts.input_file << "' for reading.\n";
    return false;
  }
  Value root;
  file >> root;
  file.close();

  data.mcell = get_node(KEY_ROOT, root, KEY_MCELL);

  // create generators
  if (data.bng_mode) {
    open_and_check_file(MODEL, bng_out, false, true);
    bng_out << GENERATED_WARNING;
    bng_gen = new BNGLGenerator(
        get_filename(data.output_files_prefix, MODEL, BNGL_EXT), bng_out, data);
  }
  python_gen = new PythonGenerator(data);

  CHECK(check_scripting(), failed);

  CHECK(generate_parameters(), failed);
  CHECK(generate_subsystem(), failed);
  std::vector<std::string> geometry_names;
  CHECK(geometry_names = generate_geometry(), failed);
  CHECK(generate_instantiation(geometry_names), failed);
  CHECK(generate_observables(opts.cellblender_viz), failed);
  CHECK(generate_model(failed), failed);
  CHECK(generate_customization(), failed);

  // delete generators
  if (data.bng_mode) {
    delete bng_gen;
    bng_gen = nullptr;
    bng_out.close();
  }
  delete python_gen;
  python_gen = nullptr;

  return !failed;
}

void MCell4Generator::check_scripting() {
  if (data.mcell.isMember(KEY_SCRIPTING)) {
    Value& scripting = get_node(data.mcell, KEY_SCRIPTING);
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
  out << GENERATED_WARNING << "\n";
  out << IMPORT_OS;
  out << MCELL_IMPORT;
  out << MODEL_PATH_SETUP << "\n";
  out << make_section_comment("model parameters");

  if (data.mcell.isMember(KEY_PARAMETER_SYSTEM)) {
    if (data.bng_mode) {
      bng_gen->generate_parameters(out);
    }
    else {
      python_gen->generate_parameters(out);
    }
  }

  out << make_section_comment("simulation setup");

  out << PARAM_ITERATIONS << " = " << data.mcell[KEY_INITIALIZATION][KEY_ITERATIONS].asString() << "\n";
  out << PARAM_TIME_STEP << " = " << data.mcell[KEY_INITIALIZATION][KEY_TIME_STEP].asString() << "\n";
  out << PARAM_DUMP << " = " << (data.debug_mode ? "True" : "False") << "\n";
  out << PARAM_EXPORT_DATA_MODEL << " = True\n";
  out << PARAM_SEED << " = 1 # may be overwritten based on commandline arguments\n";
  out << "\n";

  out.close();
}


std::string MCell4Generator::generate_species_and_mol_types(
    ostream& out, vector<SpeciesOrMolType>& species_and_mt_info) {

  // skip if there are no species/mol types
  if (!data.mcell.isMember(KEY_DEFINE_MOLECULES)) {
    return "";
  }

  if (data.bng_mode) {
    // molecule types are optional but they allow for better BNGL semantic checks

    // - we also need to generate code that sets diffusion constant,
    //   custom time/space step and other parameters to molecule types,
    // - it must be executed after bngl file is loaded, so we should not be putting it into
    //   subsystem, the new file can be called bngl_molecule_types_info.py
    // - all this extra information should be put preferably into the BNGL file
    //   but for now we would like it to be BNGL compatible and parsing comments or MCELL_ parameters is not really nice
    stringstream mt_info_out;
    bng_gen->generate_mol_types(mt_info_out);
    return mt_info_out.str();
  }
  else {
    python_gen->generate_species_and_mol_types(out, species_and_mt_info);
    return "";
  }
}


static bool rxn_uses_mcell_orientation(Value& reaction_list_item) {
  string substances =
      reaction_list_item[KEY_REACTANTS].asString() + " " +
      reaction_list_item[KEY_PRODUCTS].asString();

  size_t pos_in = substances.find(S("@") + BNG::COMPARTMENT_NAME_IN);
  size_t pos_out = substances.find(S("@") + BNG::COMPARTMENT_NAME_OUT);
  bool has_in_out_compartments = pos_in != string::npos || pos_out != string::npos;

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
      // UP is allowed when compartment class @IN or @OUT is used
      if (!(has_in_out_compartments && c == '\'')) {
        return true;
      }
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
  if (rxn_name == "") {
    return false;
  }
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


vector<IdLoc> MCell4Generator::generate_reaction_rules(ostream& out) {
  vector<IdLoc> rxn_names_w_loc;

  if (!data.mcell.isMember(KEY_DEFINE_REACTIONS)) {
    return rxn_names_w_loc;
  }

  if (data.bng_mode) {
    // put into BNG all rxn rules that,
    // 1) don't use MCell orientation or
    // 2) variable rxn rates
    // rest goes into python

    Value& define_reactions = get_node(data.mcell, KEY_DEFINE_REACTIONS);
    check_version(KEY_DEFINE_REACTIONS, define_reactions, VER_DM_2014_10_24_1638);

    bng_gen->open_reaction_rules_section();

    Value& reaction_list = get_node(define_reactions, KEY_REACTION_LIST);
    for (Value::ArrayIndex i = 0; i < reaction_list.size(); i++) {
      Value& reaction_list_item = reaction_list[i];
      check_version(KEY_MOLECULE_LIST, reaction_list_item, VER_DM_2018_01_11_1330);

      if (!rxn_uses_mcell_orientation(reaction_list_item) && !rxn_has_variable_rate(reaction_list_item)) {

        bool used_in_observables = is_rxn_used_in_observables(data.mcell, reaction_list_item[KEY_RXN_NAME].asString());
        string name = bng_gen->generate_single_reaction_rule(reaction_list_item, used_in_observables);
        rxn_names_w_loc.push_back(IdLoc(name, false));
        if (used_in_observables) {
          data.bngl_reaction_rules_used_in_observables.push_back(name);
        }
      }
      else {
        string name = python_gen->generate_single_reaction_rule(out, reaction_list_item);
        bng_gen->add_comment(
            S(IND) + "reaction '" + reaction_list_item[KEY_RXN_NAME].asString() +
            "' was generated as Python code because it contains features not supported by BNGL");
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


static void insert_compartments_in_string(const std::string& str_w_compartments, std::set<std::string>& compartments) {

  regex exp("@([0-9a-zA-Z_]+)");
  smatch res;
  string str = str_w_compartments;

  while (regex_search(str, res, exp)) {
    // not sure how else to cut off the leading '@'
    compartments.insert(string(res[0]).substr(1));
    str = res.suffix();
  }
}


void MCell4Generator::find_required_compartments(std::set<std::string>& compartments) {
  // rxns and release sites in data model may use compartments,
  // releases and observables use directly objects

  if (data.mcell.isMember(KEY_DEFINE_REACTIONS)) {
    Value& define_reactions = get_node(data.mcell, KEY_DEFINE_REACTIONS);
    Value& reaction_list = get_node(define_reactions, KEY_REACTION_LIST);
    for (Value::ArrayIndex i = 0; i < reaction_list.size(); i++) {
      Value& reaction_list_item = reaction_list[i];

      insert_compartments_in_string(reaction_list_item[KEY_REACTANTS].asString(), compartments);
      insert_compartments_in_string(reaction_list_item[KEY_PRODUCTS].asString(), compartments);
    }
  }

  if (data.mcell.isMember(KEY_RELEASE_SITES)) {
    Value& release_sites = get_node(data.mcell, KEY_RELEASE_SITES);
    Value& release_site_list = get_node(release_sites, KEY_RELEASE_SITE_LIST);
    for (Value::ArrayIndex i = 0; i < release_site_list.size(); i++) {
      Value& release_site_item = release_site_list[i];

      if (release_site_item.isMember(KEY_MOLECULE)) {
        insert_compartments_in_string(release_site_item[KEY_MOLECULE].asString(), compartments);
      }
    }
  }

  // and add all compartments used as parents & membranes
  Value& model_objects = get_node(data.mcell, KEY_MODEL_OBJECTS);
  Value& model_object_list = get_node(model_objects, KEY_MODEL_OBJECT_LIST);
  for (Value::ArrayIndex i = 0; i < model_object_list.size(); i++) {
    Value& model_object = model_object_list[i];

    const string& name = model_object[KEY_NAME].asString();
    const string& membrane_name = model_object[KEY_MEMBRANE_NAME].asString();
    const string& parent_object = model_object[KEY_PARENT_OBJECT].asString();

    if (compartments.count(name) != 0) {
      compartments.insert(membrane_name);
      compartments.insert(parent_object);
    }
  }
}


void MCell4Generator::analyze_and_generate_bngl_compartments(std::ostream& out) {

  // 1) figure out whether what compartments we need
  // not every named object must be a compartment
  find_required_compartments(data.used_compartments);

  // 2) add to BNGL file if needed
  if (data.bng_mode) {
    bng_gen->generate_compartments();
  }
}


void MCell4Generator::generate_subsystem() {
  ofstream out;
  open_and_check_file(SUBSYSTEM, out);
  out << GENERATED_WARNING << "\n";

  out << IMPORT_OS;
  out << MCELL_IMPORT;
  out << MODEL_PATH_SETUP << "\n";

  out << make_import(PARAMETERS);
  out << "\n";
  out << make_section_comment(SUBSYSTEM);

  string bngl_mol_types_initialization =
      generate_species_and_mol_types(out, data.all_species_and_mol_type_names);

  analyze_and_generate_bngl_compartments(out);

  vector<string> surface_class_names;
  python_gen->generate_surface_classes(out, surface_class_names);

  data.all_reaction_rules_names = generate_reaction_rules(out);

  out << make_section_comment("create subsystem object and add components");
  gen_ctor_call(out, SUBSYSTEM, NAME_CLASS_SUBSYSTEM, false);
  for (SpeciesOrMolType& info: data.all_species_and_mol_type_names) {
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
  for (IdLoc& r_loc: data.all_reaction_rules_names) {
    if (r_loc.in_python) {
      gen_method_call(out, SUBSYSTEM, NAME_ADD_REACTION_RULE, r_loc.name);
    }
  }

  if (data.bng_mode) {
    // each part of the bngl file is loaded separately from specific modules,
    // there is a single bngl file because later it will contain all the information needed
    // to run in it BioNetGen
    out << "\n# load subsystem information from bngl file\n";

    gen_method_call(
        out, SUBSYSTEM,
        NAME_LOAD_BNGL_MOLECULE_TYPES_AND_REACTION_RULES,
        get_abs_path(get_filename(data.output_files_prefix, MODEL, BNGL_EXT))
    );
    out << "\n" << bngl_mol_types_initialization;

    out << SET_BNGL_MOLECULE_TYPES_INFO << "(" << SUBSYSTEM << ")\n";
  }

  out.close();
}


vector<string> MCell4Generator::generate_geometry() {

  vector<string> geometry_objects;

  // TODO: check versions

  if (!data.mcell.isMember(KEY_GEOMETRICAL_OBJECTS)) {
    return geometry_objects;
  }
  Value& geometrical_objects = get_node(data.mcell, KEY_GEOMETRICAL_OBJECTS);
  if (!geometrical_objects.isMember(KEY_OBJECT_LIST)) {
    return geometry_objects;
  }

  ofstream out;
  open_and_check_file(GEOMETRY, out);
  out << GENERATED_WARNING << "\n";

  out << MCELL_IMPORT;

  python_gen->generate_geometry(out, geometry_objects);

  // NOTE: we can generate BNGL compartments from geometry here

  out.close();
  geometry_generated = true;

  return geometry_objects;
}


void MCell4Generator::generate_instantiation(const vector<string>& geometry_objects) {

  ofstream out;
  open_and_check_file(INSTANTIATION, out);
  out << GENERATED_WARNING << "\n";

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
  python_gen->generate_surface_classes_assignments(out);

  out << make_section_comment("compartments assignment");
  python_gen->generate_compartment_assignments(out);

  out << make_section_comment("create instantiation object and add components");

  gen_ctor_call(out, INSTANTIATION, NAME_CLASS_INSTANTIATION, false);
  for (const string& s: geometry_objects) {
    gen_method_call(out, INSTANTIATION, NAME_ADD_GEOMETRY_OBJECT, s);
  }
  for (const string& r: release_sites) {
    gen_method_call(out, INSTANTIATION, NAME_ADD_RELEASE_SITE, r);
  }

  out.close();
}


static string get_count_multiplier(const string& mdl_string) {

  string mdl_string_wo_comment = remove_c_comment(mdl_string);

  // check if there is a multiplier
  // handling only code generated by mcell4 for now
  string multiplier_str = "";
  size_t mult_pos = mdl_string_wo_comment.find("*");
  if (mult_pos != string::npos) {
    string substr = mdl_string_wo_comment.substr(0, mult_pos);
    try {
      stod(substr);
      multiplier_str = substr;
    }
    catch (const std::exception& e) {
      ERROR("Could not convert multiplier from " + mdl_string_wo_comment + ".");
    }
  }
  return multiplier_str;
}


void MCell4Generator::generate_counts(
    std::ostream& out, std::vector<std::string>& python_counts, bool& has_bng_observables) {

  if (!data.mcell.isMember(KEY_REACTION_DATA_OUTPUT)) {
    return;
  }

  if (data.bng_mode) {
    python_gen->generate_all_bngl_reaction_rules_used_in_observables(out);
  }

  has_bng_observables = false;

  Value& reaction_data_output = get_node(data.mcell, KEY_REACTION_DATA_OUTPUT);
  check_version(KEY_DEFINE_MOLECULES, reaction_data_output, VER_DM_2016_03_15_1800);
  string rxn_step = reaction_data_output[KEY_RXN_STEP].asString();

  Value& reaction_output_list = get_node(reaction_data_output, KEY_REACTION_OUTPUT_LIST);
  for (Value::ArrayIndex i = 0; i < reaction_output_list.size(); i++) {
    Value& reaction_output_item = reaction_output_list[i];

    string observable_name = reaction_output_item[KEY_MDL_FILE_PREFIX].asString();

    bool must_be_in_python = false;

    bool rxn_not_mol;
    bool molecules_not_species;
    bool single_term;
    string what_to_count;
    string compartment;
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
        ERROR("There is no 'COUNT' in mdl_string for output with filename " +  observable_name + ".");
      }
      else if (num_counts == 1) {
        single_term = true;
        process_single_count_term(
            data, mdl_string, rxn_not_mol, molecules_not_species,
            what_to_count, compartment, where_to_count, orientation);
      }
      else {
        must_be_in_python = true;
        single_term = false;
        what_to_count = python_gen->generate_count_terms_for_expression(out, mdl_string);
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
          observable_name + ".");
    }

    if (observable_name == "") {
      string where = where_to_count;
      if (where == "") {
        where = WORLD_FIRST_UPPER;
      }
      // using underscore instead of '.' that was used in MCell3, '.' cannot be used in BNGL
      observable_name = what_to_count + "_" + where;

      // TODO: this might need further checks
      if (observable_name.find_first_of(" ,+*/\\") != string::npos) {
        cout << "Warning: count file prefix '" + observable_name + "' is probably invalid.\n";
      }
    }

    if (data.bng_mode && !must_be_in_python &&
        bng_gen->can_express_count_with_bngl(
            single_term, rxn_not_mol, where_to_count, orientation, multiplier_str, rxn_step)) {

      if (!has_bng_observables) {
        bng_gen->open_observables_section();
      }
      bng_gen->generate_single_count(
          observable_name,
          what_to_count,
          compartment,
          molecules_not_species);
      has_bng_observables = true;
    }
    else {
      string name = create_count_name(what_to_count, compartment, where_to_count, molecules_not_species);
      python_gen->generate_single_count(
          out,
          name,
          observable_name,
          what_to_count,
          compartment,
          where_to_count, // empty for WORLD
          orientation,
          multiplier_str,
          rxn_step,
          rxn_not_mol,
          molecules_not_species,
          single_term
      );
      python_counts.push_back(name);
    }
  }
}


void MCell4Generator::generate_observables(const bool cellblender_viz) {

  ofstream out;
  open_and_check_file(OBSERVABLES, out);
  out << GENERATED_WARNING << "\n";

  out << IMPORT_OS;
  out << MCELL_IMPORT;
  out << MODEL_PATH_SETUP << "\n";

  out << make_import(PARAMETERS);
  out << make_import(SUBSYSTEM);
  if (geometry_generated) {
    out << make_import(GEOMETRY);
  }
  out << "\n";
  out << make_section_comment(OBSERVABLES);

  vector<string> viz_outputs;
  python_gen->generate_viz_outputs(out, cellblender_viz, viz_outputs);

  vector<string> counts;
  bool has_bngl_observables;
  generate_counts(out, counts, has_bngl_observables);

  out << make_section_comment("create observables object and add components");

  gen_ctor_call(out, OBSERVABLES, NAME_CLASS_OBSERVABLES, false);
  for (const string& s: viz_outputs) {
    gen_method_call(out, OBSERVABLES, NAME_ADD_VIZ_OUTPUT, s);
  }
  for (const string& r: counts) {
    gen_method_call(out, OBSERVABLES, NAME_ADD_COUNT, r);
  }
  if (has_bngl_observables) {
    out << "\n# load observables information from bngl file\n";
    gen_method_call(
        out, OBSERVABLES,
        NAME_LOAD_BNGL_OBSERVABLES,
        get_abs_path(get_filename(data.output_files_prefix, MODEL, BNGL_EXT)) + ", " +
        SUBSYSTEM + ", " +
        "'" + DEFAULT_RXN_OUTPUT_FILENAME_PREFIX + "'"
    );

    bng_gen->close_observables_section();
  }

  observables_generated = true;

  out.close();
}


void MCell4Generator::generate_config(ostream& out) {
  out << make_section_comment("configuration");

  // using values from generated parameters.py
  string params_module = get_module_name(PARAMETERS);
  gen_assign(out, MODEL, NAME_CONFIG, NAME_TIME_STEP, params_module + "." + PARAM_TIME_STEP);
  gen_assign(out, MODEL, NAME_CONFIG, NAME_SEED, params_module + "." + PARAM_SEED);
  gen_assign(out, MODEL, NAME_CONFIG, NAME_TOTAL_ITERATIONS, params_module + "." + PARAM_ITERATIONS);
  out << "\n";
  gen_assign(out, MODEL, NAME_NOTIFICATIONS, NAME_RXN_AND_SPECIES_REPORT, true);
  out << "\n";

  if (!data.mcell.isMember(KEY_INITIALIZATION)) {
    ERROR(S("Data model does not contain key ") + KEY_INITIALIZATION + ".");
  }

  Value& initialization = data.mcell[KEY_INITIALIZATION];
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
  out << GENERATED_WARNING << "\n";

  if (print_failed_marker) {
    out <<
        "ERROR: Conversion from data model failed, these generated sources are result of a "
        "best-effort conversion and will contain errors that might need to be fixed manually.\n\n";
  }

  out << BASE_MODEL_IMPORTS;
  out << get_import("importlib.util");
  out << "\n";
  out << MODEL_PATH_SETUP;
  out << "\n";
  out << MCELL_PATH_SETUP;
  out << "\n";
  out << MCELL_IMPORT;

  string parameters_module = get_module_name(PARAMETERS);
  string customization_module = get_module_name(CUSTOMIZATION);

  out <<
      "# parameters are intentionally not imported using from ... import *\n" <<
      "# because we may need to make changes to the module's variables\n" <<
      IMPORT << " " << parameters_module << "\n\n";

  out << get_resume_from_checkpoint_code(parameters_module);

  out << get_customization_import(customization_module);
  out << "\n";

  if (data.testing_mode) {
    out << CHECKPOINT_ITERATION << " = None\n\n";
  }

  out << get_argparse_w_customization_begin(parameters_module, customization_module);

  if (data.testing_mode) {
    out << get_argparse_checkpoint_iteration(parameters_module);
  }

  out << get_argparse_w_customization_end();
  out << "\n";

  out << get_import(get_module_name(SUBSYSTEM));
  out << get_import(get_module_name(INSTANTIATION));
  if (observables_generated) {
    out << get_import(get_module_name(OBSERVABLES));
  }
  out << "\n";

  out << "# create main model object\n";
  gen_ctor_call(out, MODEL, NAME_CLASS_MODEL, false);
  out << "\n";

  generate_config(out);

  if (data.testing_mode) {
    out << make_section_comment("testing-specific configuration");
    gen_assign(out, MODEL, NAME_CONFIG, NAME_REACTION_CLASS_CLEANUP_PERIODICITY, TESTING_RXN_CLASS_CLEANUP_PERIODICITY);
    gen_assign(out, MODEL, NAME_CONFIG, NAME_SPECIES_CLEANUP_PERIODICITY, TESTING_SPECIES_CLEANUP_PERIODICITY);
    out << "\n";
  }

  out << get_user_defined_configuration(customization_module);

  out << "\n";

  out << make_section_comment("add components");
  gen_method_call(out, MODEL, NAME_ADD_SUBSYSTEM, get_module_name(SUBSYSTEM) + "." + SUBSYSTEM);
  gen_method_call(out, MODEL, NAME_ADD_INSTANTIATION, get_module_name(INSTANTIATION) + "." + INSTANTIATION);
  if (observables_generated) {
    gen_method_call(out, MODEL, NAME_ADD_OBSERVABLES, get_module_name(OBSERVABLES) + "." + OBSERVABLES);
  }
  out << "\n";

  out << make_section_comment("initialization and execution");

  out <<
      "if " << customization_module << " and '" << CUSTOM_INIT_AND_RUN << "' in dir(" << customization_module << "):\n" <<
      IND4 << customization_module << "." << CUSTOM_INIT_AND_RUN << "(" << MODEL << ")\n" <<
      "else:\n"
  ;

  out << IND4;
  gen_method_call(out, MODEL, NAME_INITIALIZE);
  out << "\n";

  if (!data.checkpoint_iterations.empty()) {
    out << make_section_comment("checkpoint iterations");
    for (int it: data.checkpoint_iterations) {
      out << IND4;
      gen_method_call(out, MODEL, NAME_SCHEDULE_CHECKPOINT, to_string(it));
    }
    out << "\n";
  }

  out << IND4 << "if " << parameters_module << "." << PARAM_DUMP << ":\n";
  out << IND8;
  gen_method_call(out, MODEL, NAME_DUMP_INTERNAL_STATE);
  out << "\n";

  // method export_data_model uses target directory from viz_outputs
  out << IND4 << "if " << parameters_module << "." << PARAM_EXPORT_DATA_MODEL << " and " << MODEL << "." << NAME_VIZ_OUTPUTS << ":\n";
  out << IND8;
  gen_method_call(out, MODEL, NAME_EXPORT_DATA_MODEL);
  out << "\n";

  out << IND4;
  gen_method_call(out, MODEL, NAME_RUN_ITERATIONS, parameters_module + "." + PARAM_ITERATIONS);
  out << IND4;
  gen_method_call(out, MODEL, NAME_END_SIMULATION);
}


void MCell4Generator::generate_customization() {
  // check if file exists, do not overwrite
  string cust_filename = get_filename(data.output_files_prefix, CUSTOMIZATION, PY_EXT);
  ifstream tmp;
  tmp.open(cust_filename);
  if (tmp.is_open()) {
    // TODO: test
    cout << "Message: custom file " + cust_filename + " already exists, keeping it as it is.\n";
    tmp.close();
    return;
  }

  ofstream out;
  open_and_check_file(CUSTOMIZATION, out);

  out <<
      "# This file hooks to override default MCell4 model\n"
      "# code behavior for models generated from CellBlender\n";

  out << MCELL_IMPORT;
  out << IMPORT << " " << get_module_name(PARAMETERS) << "\n\n";

  out << TEMPLATE_CUSTOM_ARGPARSE_AND_PARAMETERS << "\n";
  out << TEMPLATE_CUSTOM_CONFIG << "\n";
  out << TEMPLATE_CUSTOM_INIT_AND_RUN << "\n";
}


} // namespace MCell
