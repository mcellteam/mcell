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
    bng_gen->generate_units_information_header();
  }
  python_gen = new PythonGenerator(data);


  CHECK(generate_customization(), failed);

  CHECK(generate_shared(), failed);
  CHECK(generate_parameters(), failed);
  CHECK(generate_subsystem(), failed);
  std::vector<std::string> geometry_names;
  CHECK(geometry_names = generate_geometry(), failed);
  CHECK(generate_instantiation(geometry_names), failed);
  CHECK(generate_observables(), failed);
  CHECK(generate_model(failed), failed);
  CHECK(generate_customization_template(), failed);

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


static std::string get_file_base_name(const std::string& path) {
  size_t pos = path.find_last_of("/\\");
  if (pos != string::npos) {
    return path.substr(pos + 1);
  }
  else {
    return path;
  }
}


void MCell4Generator::generate_customization() {
  // first check unsupported MCell3 scripting
  Value& scripting = get_node(data.mcell, KEY_SCRIPTING);
  if (scripting.isMember(KEY_SCRIPTING_LIST)) {
    Value& scripting_list = get_node(scripting, KEY_SCRIPTING_LIST);
    if (!scripting_list.empty()) {
      ERROR("Data model contains MCell3 scripting. To convert this model: 1) disable MCell3 mode, 2) generate MDL, "
          "3) convert MDL to data model and 4) use this converter. "
          "Conversion will now continue, but it will be missing the scripted code."
      );
    }
  }

  // copy files from MCell4 scripting
  if (scripting.isMember(KEY_MCELL4_SCRIPTING_LIST)) {
    Value& mcell4_scripting_list = get_node(scripting, KEY_MCELL4_SCRIPTING_LIST);
    for (Value::ArrayIndex i = 0; i < mcell4_scripting_list.size(); i++) {
      Value& item = mcell4_scripting_list[i];
      bool internal;
      string internal_external = get_node(item, KEY_INTERNAL_EXTERNAL).asString();
      if (internal_external == VALUE_INTERNAL) {
        internal = true;
      }
      else if (internal_external == VALUE_EXTERNAL) {
        internal = false;
      }
      else {
        ERROR(S(KEY_INTERNAL_EXTERNAL) + " may be either " + VALUE_INTERNAL + " or " + VALUE_EXTERNAL +
            ", not " + internal_external + ".");
      }

      if (internal) {
        // create a file from internal script
        string fname = get_node(item, KEY_INTERNAL_FILE_NAME).asString();
        Value& script_texts = get_node(scripting, KEY_SCRIPT_TEXTS);
        if (!script_texts.isMember(fname)) {
          ERROR("Internal script " + fname + " was not found.");
        }
        // expecting the the file can be represented as text (it must be since it is in a JSON format)
        string code = get_node(script_texts, fname).asString();
        ofstream fout;
        fout.open(fname);
        if (!fout.is_open()) {
          ERROR("Could not open file " + fname + " for writing for internal script export.");
        }

        cout << "Creating " + fname + " from internal script file.\n";
        fout.write(code.c_str(), code.size());
        fout.close();
      }
      else {
        // copy external file
        string fname = get_node(item, KEY_EXTERNAL_FILE_NAME).asString();

        ifstream fin;
        fin.open(fname);
        if (!fin.is_open()) {
          ERROR("Could not open file " + fname + " for reading for external script export.");
        }

        // get base path
        string exported_fname = get_file_base_name(fname);
        cout << "exported_fname: " << exported_fname << "\n";
        ofstream fout;
        fout.open(exported_fname);
        if (!fout.is_open()) {
          fin.close();
          ERROR("Could not open file " + fname + " for writing for external script export.");
        }

        // copy data
        cout << "Copying " + fname + " to " + exported_fname + ".\n";
        fout << fin.rdbuf();
        fin.close();
        fout.close();
      }
    }
  }
}


void MCell4Generator::generate_shared() {
  ofstream out;
  open_and_check_file_w_prefix("", SHARED, out);
  out << GENERATED_WARNING << "\n";
  out <<
      "# This is an auxiliary module containing only a dictionary for parameter overrides.\n" <<
      "# It has to be a separate module so that it can be shared easily among all modules.\n";
  out << PARAMETER_OVERRIDES << " = {}\n";
  out.close();
}


void MCell4Generator::generate_simulation_setup_parameter(
    std::ostream& out, const string& name, const string& value) {
  string ind;
  if (!data.not_overridable_python_params) {
    // allow params to be overridden
    out << "if " << NOT_DEFINED << "('" << name << "'):\n";
    ind = IND4;
  }
  out << ind << name << " = " << value << "\n\n";
}


void MCell4Generator::generate_parameters() {
  ofstream out;
  open_and_check_file(PARAMETERS, out);
  out << GENERATED_WARNING << "\n";
  out << IMPORT_SYS_OS;
  out << IMPORT_MATH;
  out << IMPORT_SHARED;
  out << IMPORT_MCELL_AS_M;
  out << MODEL_PATH_SETUP << "\n";
  out << make_section_comment("model parameters");

  out <<
    "# declare all items from parameter_overrides as variables\n" <<
    "for parameter_name, value in " << SHARED << "." << PARAMETER_OVERRIDES << ".items():\n" <<
    "    setattr(sys.modules[__name__], parameter_name, value)\n" <<
    "\n" <<
    "# auxiliary function used to determine whether a parameter was defined\n" <<
    "def " << NOT_DEFINED << "(parameter_name):\n" <<
    "    return parameter_name not in globals()\n\n";

  if (data.mcell.isMember(KEY_PARAMETER_SYSTEM)) {
    if (data.bng_mode) {
      bng_gen->generate_parameters(out);
    }
    else {
      python_gen->generate_parameters(out);
    }
  }

  out << make_section_comment("simulation setup");

  generate_simulation_setup_parameter(out, PARAM_ITERATIONS, data.mcell[KEY_INITIALIZATION][KEY_ITERATIONS].asString());
  generate_simulation_setup_parameter(out, PARAM_TIME_STEP, data.mcell[KEY_INITIALIZATION][KEY_TIME_STEP].asString());
  generate_simulation_setup_parameter(out, PARAM_DUMP, (data.debug_mode ? "True" : "False"));
  generate_simulation_setup_parameter(out, PARAM_EXPORT_DATA_MODEL, "True");
  generate_simulation_setup_parameter(out, PARAM_SEED, "1");
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


static bool rxn_uses_surf_class(
    Value& reaction_list_item,
    const std::vector<std::string>& surf_class_names) {

  vector<string> substances;
  vector<string> orientations;
  Value& reactants = reaction_list_item[KEY_REACTANTS];
  parse_rxn_rule_side(reactants, substances, orientations);

  for (const string& reac: substances) {
    if (find(surf_class_names.begin(), surf_class_names.end(), reac) != surf_class_names.end()) {
      return true;
    }
  }
  return false;
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


vector<IdLoc> MCell4Generator::generate_reaction_rules(ostream& out, const std::vector<std::string>& surf_class_names) {
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

      if (!rxn_uses_mcell_orientation(reaction_list_item) &&
          !rxn_has_variable_rate(reaction_list_item) &&
          !rxn_uses_surf_class(reaction_list_item, surf_class_names)) {

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
  out << IMPORT_SHARED;
  out << IMPORT_MCELL_AS_M;
  out << make_import(PARAMETERS);
  out << "\n";
  out << make_section_comment(SUBSYSTEM);
  out << MODEL_PATH_SETUP << "\n";

  string bngl_mol_types_initialization =
      generate_species_and_mol_types(out, data.all_species_and_mol_type_names);

  analyze_and_generate_bngl_compartments(out);

  vector<string> surface_class_names;
  python_gen->generate_surface_classes(out, surface_class_names);

  data.all_reaction_rules_names = generate_reaction_rules(out, surface_class_names);

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
        get_abs_path(get_filename(data.output_files_prefix, MODEL, BNGL_EXT)),
        S(SHARED) + "." + PARAMETER_OVERRIDES
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

  out << IMPORT_MCELL_AS_M;

  python_gen->generate_geometry(out, geometry_objects);

  // NOTE: we can generate BNGL compartments from geometry here

  out.close();
  geometry_generated = true;

  return geometry_objects;
}


void MCell4Generator::generate_release_sites(std::ostream& out, std::vector<std::string>& release_site_names) {

  if (!data.mcell.isMember(KEY_RELEASE_SITES)) {
    return;
  }

  Value& release_sites = data.mcell[KEY_RELEASE_SITES];
  check_version(KEY_RELEASE_SITES, release_sites, VER_DM_2014_10_24_1638);
  Value& release_site_list = release_sites[KEY_RELEASE_SITE_LIST];
  if (release_site_list.empty()) {
    return;
  }

  if (data.bng_mode) {
    bool can_express_all_releases_with_bngl = true;

    for (Value::ArrayIndex i = 0; i < release_site_list.size(); i++) {
      // simulation result differs based on the order of releases so we must either generate all
      // releases into BNGL or none
      Value& release_site_item = release_site_list[i];
      if (!bng_gen->can_express_release_with_bngl(release_site_item)) {
        can_express_all_releases_with_bngl = false;
        break;
      }
    }
    if (can_express_all_releases_with_bngl) {
      bng_gen->open_seed_species_section();

      for (Value::ArrayIndex i = 0; i < release_site_list.size(); i++) {
        // simulation result differs based on the order of releases so we must either generate all
        // releases into BNGL or none
        const Value& release_site_item = release_site_list[i];
        check_not_empty(release_site_item, KEY_QUANTITY, "Release site");
        bng_gen->generate_single_release_site(
            release_site_item[KEY_MOLECULE].asString(),
            release_site_item[KEY_QUANTITY].asString(),
            get_description(release_site_item));
      }

      bng_gen->close_seed_species_section();

      // all done
      return;
    }
  }

  // python variant - used either without bng mode or as a fallback for bng mode
  python_gen->generate_release_sites(out, release_site_names);
}


void MCell4Generator::generate_instantiation(const vector<string>& geometry_objects) {

  ofstream out;
  open_and_check_file(INSTANTIATION, out);
  out << GENERATED_WARNING << "\n";

  out << IMPORT_OS;
  out << IMPORT_SHARED;
  out << IMPORT_MCELL_AS_M;
  out << make_import(PARAMETERS);
  out << make_import(SUBSYSTEM);
  if (geometry_generated) {
    out << make_import(GEOMETRY);
  }
  out << MODEL_PATH_SETUP << "\n";
  out << "\n";
  out << make_section_comment(INSTANTIATION);
  out << make_section_comment("release sites");

  out << make_section_comment("surface classes assignment");
  python_gen->generate_surface_classes_assignments(out);

  out << make_section_comment("compartments assignment");
  python_gen->generate_compartment_assignments(out);

  vector<string> release_sites;
  generate_release_sites(out, release_sites);

  out << make_section_comment("create instantiation object and add components");

  gen_ctor_call(out, INSTANTIATION, NAME_CLASS_INSTANTIATION, false);
  for (const string& s: geometry_objects) {
    gen_method_call(out, INSTANTIATION, NAME_ADD_GEOMETRY_OBJECT, s);
  }
  for (const string& r: release_sites) {
    gen_method_call(out, INSTANTIATION, NAME_ADD_RELEASE_SITE, r);
  }

  if (data.bng_mode) {
    out << "\n# load seed species information from bngl file\n";

    gen_method_call(
        out, INSTANTIATION,
        NAME_LOAD_BNGL_COMPARTMENTS_AND_SEED_SPECIES,
        get_abs_path(get_filename(data.output_files_prefix, MODEL, BNGL_EXT)),
        data.has_default_compartment_object ? BNG::DEFAULT_COMPARTMENT_NAME : "None",
        S(SHARED) + "." + PARAMETER_OVERRIDES
    );
    out << "\n";
  }

  out.close();
}


static size_t find_closing_bracket_pos(const string& s, const size_t start_pos) {
  if (start_pos == string::npos) {
    return string::npos;
  }

  assert(s[start_pos] == '[');

  size_t pos = start_pos + 1;
  int bracket_count = 1;
  while (pos < s.size() && bracket_count != 0) {
    switch (s[pos]) {
      case '[':
        bracket_count++;
        break;
      case ']':
        bracket_count--;
        break;
      default:
        // skip
        break;
    }
    pos++;
  }
  if (bracket_count != 0) {
    return string::npos;
  }
  else {
    return pos - 1;
  }
}


// also checks that the MDLString is correctly formed
static string check_mdlstring_count_and_get_mult_or_div(const string& mdl_string) {

  const string format_msg =
      ". MDL string supported by MCell4 must be in format 'COUNT[pattern,region] *|/ const_expr' or "
      "{(} COUNT[pattern,region] { +|- {(} COUNT[pattern,region] {)} } {)} *|/ const_expr' "
      "({x} means 0..n of repetitions of x, and x|y means x or y).";

  string mdl_string_wo_comment = remove_c_comment(mdl_string);
  string str = remove_whitespace(mdl_string_wo_comment);

  int num_parens = 0;
  size_t pos = 0;
  const string COUNT = "COUNT";

  // check until we process all COUNTs with all parentheses
  while (pos < str.size() && !(str.find(COUNT, pos) == string::npos && num_parens == 0)) {

    switch (str[pos]) {
      case 'C': {
          // this must be COUNT[...]
          size_t count_pos = str.find(COUNT, pos);
          size_t opening_bracket_pos = str.find('[', pos);
          // brackets are allowed inside e,g, such as here: COUNT[sm1,Scene.Cube[ALL]]
          size_t closing_bracket_pos = find_closing_bracket_pos(str, opening_bracket_pos);
          if (pos != count_pos || opening_bracket_pos == string::npos || closing_bracket_pos == string::npos) {
            ERROR("Could not parse COUNT specifier in " + mdl_string_wo_comment + format_msg);
          }
          pos = closing_bracket_pos + 1;
        }
        break;

      case '(':
        num_parens++;
        pos++;
        break;

      case ')':
        num_parens--;
        pos++;
        break;

      case '+':
      case '-':
        // allowed
        pos++;
        break;

      default:
        ERROR("Unexpected character '" + str[pos] + "' in MDL string " +
            mdl_string_wo_comment + format_msg);
        break;
    }
  }

  if (num_parens != 0) {
    ERROR("Not matching parentheses in MDL string " +
        mdl_string_wo_comment + format_msg);
  }

  if (pos == str.size()) {
    // all checked, nothing more to process
    return "";
  }

  if (pos < str.size() && str.find(COUNT, pos) != string::npos) {
    ERROR("Unexpected COUNT specifier in MDL string " +
        mdl_string_wo_comment + format_msg);
  }


  // what optionally follows must be * or /
  string mul_or_div_expr_str = str.substr(pos);
  char first_expr_char = mul_or_div_expr_str[0];
  if (first_expr_char != '*' && first_expr_char != '/') {
    ERROR("The only allowed operator following COUNT expression(s) is '*' or '/', error for '" + first_expr_char + "' in " +
        mdl_string_wo_comment + format_msg);
  }

  return mul_or_div_expr_str;
}


// sets is_gdat, may return empty string meaning that the name must be determined automatically
static string get_count_file_name(Value& reaction_output_item, bool& is_gdat) {
  if (reaction_output_item.isMember(KEY_OUTPUT_FILE_OVERRIDE)) {

    string res = reaction_output_item[KEY_OUTPUT_FILE_OVERRIDE].asString();
    string gdat = "gdat";
    is_gdat = res.size() > gdat.size() && res.substr(res.size() - gdat.size()) == gdat;
    return res;
  }
  else {
    // default
    string mdl_file_prefix = reaction_output_item[KEY_MDL_FILE_PREFIX].asString();

    is_gdat = false;

    if (mdl_file_prefix == "") {
      // not set
      return "";
    }
    else {
      return mdl_file_prefix + ".dat";
    }
  }
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
    bool is_gdat = false;
    string file_name = get_count_file_name(reaction_output_item, is_gdat);

    bool must_be_in_python = false;

    bool rxn_not_mol;
    bool molecules_not_species;
    bool single_term;
    string what_to_count;
    string where_to_count; // empty for WORLD
    string orientation;

    string count_location = reaction_output_item[KEY_COUNT_LOCATION].asString();

    string rxn_or_mol = reaction_output_item[KEY_RXN_OR_MOL].asString();
    string mul_div_str = "";
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
            what_to_count, where_to_count, orientation);
      }
      else {
        must_be_in_python = true;
        single_term = false;
      }

      mul_div_str = check_mdlstring_count_and_get_mult_or_div(mdl_string);
      assert(mul_div_str.empty() || mul_div_str[0] == '*' || mul_div_str[0] == '/');
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
      // this is a case where mdl_string is empty

      string where = where_to_count;
      if (where == "") {
        where = WORLD_FIRST_UPPER;
      }
      // using underscore instead of '.' that was used in MCell3, '.' cannot be used in BNGL
      observable_name = what_to_count + "_" + where;
      if (mul_div_str != "") {
        observable_name += '_' + fix_id(mul_div_str);
      }

      // TODO: this might need further checks
      if (observable_name.find_first_of(" ,+*/\\") != string::npos) {
        cout << "Warning: count file prefix '" + observable_name + "' is probably invalid.\n";
      }

      if (file_name == "") {
        file_name = observable_name + ".dat";
      }
    }

    if (data.bng_mode && !must_be_in_python &&
        bng_gen->can_express_count_with_bngl(
            single_term, rxn_not_mol, where_to_count, orientation, mul_div_str, rxn_step)) {

      if (!has_bng_observables) {
        bng_gen->open_observables_section();
      }
      bng_gen->generate_single_count(
          observable_name,
          what_to_count,
          get_description(reaction_output_item),
          molecules_not_species);
      has_bng_observables = true;

      if (is_gdat) {
        if (data.bng_observables_output_gdat_file == "") {
          data.bng_observables_output_gdat_file = file_name;
        }
        else if (file_name != data.bng_observables_output_gdat_file) {
          ERROR("Cannot use multiple GDAT output files with BNGL mode, error for " + file_name);
        }
      }
    }
    else {
      string name = fix_id(COUNT_PREFIX + observable_name);

      // prepare count terms
      string mdl_string = reaction_output_item[KEY_MDL_STRING].asString();
      string count_term_name =
          python_gen->generate_count_terms_for_expression(
              out, mdl_string, what_to_count, where_to_count, orientation, rxn_not_mol);
      where_to_count = "";

      gen_description(out, reaction_output_item);

      python_gen->generate_single_count(
          out,
          name,
          observable_name,
          file_name,
          count_term_name,
          mul_div_str,
          rxn_step
      );
      python_counts.push_back(name);
    }
  }
}


void MCell4Generator::generate_observables() {

  ofstream out;
  open_and_check_file(OBSERVABLES, out);
  out << GENERATED_WARNING << "\n";

  out << IMPORT_OS;
  out << IMPORT_SHARED;
  out << IMPORT_MCELL_AS_M;
  out << make_import(PARAMETERS);
  out << make_import(SUBSYSTEM);
  if (geometry_generated) {
    out << make_import(GEOMETRY);
  }
  out << MODEL_PATH_SETUP << "\n";
  out << "\n";
  out << make_section_comment(OBSERVABLES);

  bool use_cellblender_output;
  if (data.mcell[KEY_INITIALIZATION].isMember(KEY_EXPORT_ALL_ASCII)) {
    use_cellblender_output = !data.mcell[KEY_INITIALIZATION][KEY_EXPORT_ALL_ASCII].asBool();
  }
  else {
    use_cellblender_output = false;
  }

  vector<string> viz_outputs;
  python_gen->generate_viz_outputs(out, use_cellblender_output, viz_outputs);

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
        "'" + DEFAULT_RXN_OUTPUT_FILENAME_PREFIX + data.bng_observables_output_gdat_file + "'",
        S(SHARED) + "." + PARAMETER_OVERRIDES
    );

    bng_gen->close_observables_section();
  }

  observables_generated = true;

  out.close();
}

// returns true for "ON and false for "OFF, fails when a different value is used
static std::string convert_warning_level(const std::string& value) {
  if (value == VALUE_ERROR) {
    return S(MDOT) + NAME_ENUM_WARNING_LEVEL + "." + NAME_EV_ERROR;
  }
  else if (value == VALUE_WARNING) {
    return S(MDOT) + NAME_ENUM_WARNING_LEVEL + "." + NAME_EV_WARNING;
  }
  else if (value == VALUE_IGNORED) {
    return S(MDOT) + NAME_ENUM_WARNING_LEVEL + "." + NAME_EV_IGNORE;
  }
  else {
    ERROR("Invalid value " + value + ", expected only WARNING or IGNORED.\n");
  }
}


void MCell4Generator::generate_config(ostream& out) {
  out << make_section_comment("configuration");

  // using values from generated parameters.py where applicable
  gen_assign(out, MODEL, NAME_CONFIG, NAME_TIME_STEP, PARAM_TIME_STEP);
  if (data.mcell.isMember(KEY_USE_BNG_UNITS) && data.mcell[KEY_USE_BNG_UNITS].asBool()) {
    gen_assign(out, MODEL, NAME_CONFIG, NAME_USE_BNG_UNITS, true);
  }
  gen_assign(out, MODEL, NAME_CONFIG, NAME_SEED, PARAM_SEED);
  gen_assign(out, MODEL, NAME_CONFIG, NAME_TOTAL_ITERATIONS, PARAM_ITERATIONS);
  out << "\n";

  Value& initialization = data.mcell[KEY_INITIALIZATION];

  if (initialization.isMember(KEY_WARNINGS)) {
    Value& warnings = initialization[KEY_WARNINGS];

    if (warnings.isMember(KEY_HIGH_REACTION_PROBABILITY)) {
      gen_assign(out, MODEL, NAME_WARNINGS, NAME_HIGH_REACTION_PROBABILITY,
          convert_warning_level(warnings[KEY_HIGH_REACTION_PROBABILITY].asString())
      );
    }
    if (warnings.isMember(KEY_MOLECULE_PLACEMENT_FAILURE)) {
      gen_assign(out, MODEL, NAME_WARNINGS, NAME_MOLECULE_PLACEMENT_FAILURE,
          convert_warning_level(warnings[KEY_HIGH_REACTION_PROBABILITY].asString())
      );
    }
  }

  if (initialization.isMember(KEY_NOTIFICATIONS)) {
    Value& notifications = initialization[KEY_NOTIFICATIONS];

    if (notifications.isMember(KEY_SPECIES_REACTIONS_REPORT)) {
      gen_assign(out, MODEL, NAME_NOTIFICATIONS, NAME_RXN_AND_SPECIES_REPORT,
          notifications[KEY_SPECIES_REACTIONS_REPORT].asBool()
      );
    }

    if (notifications.isMember(KEY_VARYING_PROBABILITY_REPORT)) {
      gen_assign(out, MODEL, NAME_NOTIFICATIONS, NAME_RXN_PROBABILITY_CHANGED,
          notifications[KEY_VARYING_PROBABILITY_REPORT].asBool()
      );
    }
  }

  out << "\n";

  if (!data.mcell.isMember(KEY_INITIALIZATION)) {
    ERROR(S("Data model does not contain key ") + KEY_INITIALIZATION + ".");
  }

  Value& partitions = initialization[KEY_PARTITIONS];

  // choose the largest value for partition size and the smallest step
  double x_start = stod(partitions[KEY_X_START].asString());
  double x_end = stod(partitions[KEY_X_END].asString());
  double x_step = stod(partitions[KEY_X_STEP].asString());
  double y_start = stod(partitions[KEY_Y_START].asString());
  double y_end = stod(partitions[KEY_Y_END].asString());
  double y_step = stod(partitions[KEY_Y_STEP].asString());
  double z_start = stod(partitions[KEY_Z_START].asString());
  double z_end = stod(partitions[KEY_Z_END].asString());
  double z_step = stod(partitions[KEY_Z_STEP].asString());

  double partition_dimension;
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


  double partition_step;
  if (!cmp_eq(x_step, y_step) || !cmp_eq(y_step, z_step)) {
    cout << "Message: Individual partition step sizes are different, changing the step to be the smallest of them.\n";

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
    gen_assign(out, MODEL, NAME_CONFIG, NAME_INTERACTION_RADIUS, radius);
  }

  string surf_grid_density = initialization[KEY_SURFACE_GRID_DENSITY].asString();
  if (surf_grid_density != "" && surf_grid_density != "10000") {
    gen_assign(out, MODEL, NAME_CONFIG, NAME_SURFACE_GRID_DENSITY, surf_grid_density);
  }

  string vacancy_search_distance = initialization[KEY_VACANCY_SEARCH_DISTANCE].asString();
  if (vacancy_search_distance != "" && vacancy_search_distance != "10") {
    gen_assign(out, MODEL, NAME_CONFIG, NAME_VACANCY_SEARCH_DISTANCE, vacancy_search_distance);
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

  out << IMPORT_SYS_OS;
  out << get_import("importlib.util");
  out << "\n";
  out << MODEL_PATH_SETUP;
  out << "\n";
  out << MCELL_PATH_SETUP;
  out << "\n";
  out << IMPORT_MCELL_AS_M;

  string customization_module = CUSTOMIZATION;
  string shared_module = SHARED;

  out << "\n" << make_section_comment("customization and argument processing");
  out << "# this module is used to hold any overrides of parameter values\n";
  out << IMPORT_SHARED;

  out << "# import the customization.py module if it exists\n";
  out << get_customization_import(customization_module);
  out << "\n";

  if (data.testing_mode) {
    out << CHECKPOINT_ITERATION << " = None\n\n";
  }

  out << "# process command-line arguments\n";
  out << get_argparse_w_customization_begin(customization_module);
  if (data.testing_mode) {
    out << get_argparse_checkpoint_iteration();
  }

  out << get_argparse_w_customization_end();
  out << "\n";

  out << S("\n# the module parameters uses ") + SHARED + "." + PARAMETER_OVERRIDES + " to override parameter values\n";
  out << make_import(PARAMETERS);
  out << "\n\n";

  out << get_resume_from_checkpoint_code();


  out << "\n" << make_section_comment("model creation and simulation");

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

  out << make_section_comment("add components");
  out << get_import(get_module_name(SUBSYSTEM));
  out << get_import(get_module_name(INSTANTIATION));
  if (observables_generated) {
    out << get_import(get_module_name(OBSERVABLES));
  }
  out << "\n";

  gen_method_call(out, MODEL, NAME_ADD_SUBSYSTEM, get_module_name(SUBSYSTEM) + "." + SUBSYSTEM);
  gen_method_call(out, MODEL, NAME_ADD_INSTANTIATION, get_module_name(INSTANTIATION) + "." + INSTANTIATION);
  if (observables_generated) {
    gen_method_call(out, MODEL, NAME_ADD_OBSERVABLES, get_module_name(OBSERVABLES) + "." + OBSERVABLES);
  }
  out << "\n";

  out << get_user_defined_configuration(customization_module);

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

  out << IND4 << "if " << PARAM_DUMP << ":\n";
  out << IND8;
  gen_method_call(out, MODEL, NAME_DUMP_INTERNAL_STATE);
  out << "\n";

  // method export_data_model uses target directory from viz_outputs
  out << IND4 << "if " << PARAM_EXPORT_DATA_MODEL << " and " << MODEL << "." << NAME_VIZ_OUTPUTS << ":\n";
  out << IND8;
  gen_method_call(out, MODEL, NAME_EXPORT_DATA_MODEL);
  out << "\n";

  out << IND4;
  gen_method_call(out, MODEL, NAME_RUN_ITERATIONS, PARAM_ITERATIONS);
  out << IND4;
  gen_method_call(out, MODEL, NAME_END_SIMULATION);
}


void MCell4Generator::generate_customization_template() {
  // check if file exists, do not overwrite
  string cust_filename = S(CUSTOMIZATION) + PY_EXT;
  ifstream tmp;
  tmp.open(cust_filename);
  if (tmp.is_open()) {
    cout << "Message: custom file " + cust_filename + " already exists, keeping it as it is.\n";
    tmp.close();
    return;
  }

  ofstream out;
  open_and_check_file_w_prefix("", CUSTOMIZATION, out);

  out <<
      "# This file contains hooks to override default MCell4 model\n"
      "# code behavior for models generated from CellBlender\n";

  out << IMPORT_SYS_OS;
  out << IMPORT_SHARED;
  out << IMPORT_MCELL_AS_M;

  out << TEMPLATE_CUSTOM_ARGPARSE_AND_PARAMETERS << "\n";
  out << TEMPLATE_CUSTOM_CONFIG << "\n";
  out << get_template_custom_init_and_run(get_module_name(PARAMETERS)) << "\n";
}


} // namespace MCell
