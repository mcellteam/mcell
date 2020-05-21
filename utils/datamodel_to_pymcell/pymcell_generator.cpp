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

#include <datamodel_to_pymcell/pymcell_generator.h>
#include <exception>
#include <iostream>
#include <string>
#include <cassert>

#include "../../libmcell/generated/gen_names.h"
#include "../../include/datamodel_defines.h"

using namespace std;
using namespace MCell::API;

namespace MCell {

using Json::Value;

typedef std::invalid_argument ConversionError;


const char* const PARAMETERS_SUFFIX = "parameters";
const char* const SUBSYSTEM_SUFFIX = "subsystem";
const char* const GEOMETRY_SUFFIX = "geometry";
const char* const INSTANTIATION_SUFFIX = "instantiation";
const char* const OBSERVABLES_SUFFIX = "observables";
const char* const MODEL_SUFFIX = "model";

const char* const PY_EXT = ".py";

const char* const MDOT = "m.";
const char* const MCELL_IMPORT = "import mcell as m\n\n";

const char* const IND = "    ";
const char* const BLOCK_BEGIN1 = "# ---- ";
const char* const BLOCK_BEGIN2 = " ----\n";
const char* const BLOCK_END1 = "# ^^^^ ";
const char* const BLOCK_END2 = " ^^^^\n\n";
const char* const CTOR_END = ")\n\n";

const char* const PARAM_SEED = "SEED";
const char* const PARAM_ITERATIONS = "ITERATIONS";
const char* const PARAM_TIME_STEP = "TIME_STEP";
const char* const PARAM_DUMP = "DUMP";

const char* const UNNAMED_REACTION_RULE_PREFIX = "unnamed_reaction_rule_";

// using exceptions to

#define CHECK(stmt, failed) \
  do { \
    try { \
      (stmt); \
    } \
    catch (exception& e) { \
      cerr << e.what() << "\n"; \
      cerr << "Exception caught in '" << __FUNCTION__ << "' after conversion error.\n"; \
      failed = true; \
    } \
  } while (0)

#define CHECK_PROPERTY(cond) \
  do { \
    if(!(cond)) { \
      throw ConversionError(S("Expected '") + #cond + "' is false. (" + __FUNCTION__ + " - " + __FILE__ + ":" + to_string(__LINE__) + ")"); \
    } \
  } while (0)

#define ERROR(msg) throw ConversionError(msg);


// auxiliary method to simply convert to std::string for when concatenating string
static string S(const char* s) {
  return string(s);
}


// throws exception when the member is member is there
static Value& get_node(const string parent_name, Value& parent, const string name) {
  if (!parent.isMember(name)) {
    throw ConversionError("Node '" + parent_name + "' does not contain expected node '" + name + "'.");
  }
  return parent[name];
}


// used when we know that the member is there
static Value& get_node(Value& parent, const string name) {
  assert(parent.isMember(name));
  return parent[name];
}


static void print_comma(ostream& out, Value::ArrayIndex index, Value& array) {
  if (index + 1 != array.size()) {
    out << ", ";
  }
}


template<typename T>
void print_comma(ostream& out, size_t index, const vector<T>& array) {
  if (index + 1 != array.size()) {
    out << ", ";
  }
}


static string make_section_comment(const string text) {
  return BLOCK_BEGIN1 + text + BLOCK_BEGIN2 + "\n";
}


static string make_start_block_comment(const string text) {
  return BLOCK_BEGIN1 + text + BLOCK_BEGIN2;
}


static string make_end_block_comment(const string text) {
  return BLOCK_END1 + text + BLOCK_END2 + "\n";
}


static void check_version(const string node_name, Json::Value& node, const char* const version) {
  if (node[KEY_DATA_MODEL_VERSION].asString() != version) {
    throw ConversionError(
        "Error: version for " + node_name + " is " + node[KEY_DATA_MODEL_VERSION].asString() +
        ", expected" + version);
  }
}


template<typename T>
void gen_param(ofstream& out, std::string name, T value, bool comma) {
  out << IND << name << " = " << value << (comma?",":"") << "\n";
}


template<>
void gen_param(ofstream& out, std::string name, Json::Value& value, bool comma) {
  out << IND << name << " = \"" << value.asString() << "\"" << (comma?",":"") << "\n";
}


template<>
void gen_param(ofstream& out, std::string name, std::string value, bool comma) {
  out << IND << name << " = \"" << value << "\"" << (comma?",":"") << "\n";
}


static void gen_rxn_substance_inst(ofstream& out, Json::Value& substances_node) {
  string str = substances_node.asString();

  vector<string> substances;

  // finite automata to parse the reaction side string, e.g. "a + b"
  enum state_t {
    START,
    ID,
    AFTER_ID,
    AFTER_PLUS
  };

  state_t state = START;
  string current_id;
  for (size_t i = 0; i < str.size(); i++) {
    char c = str[i];
    switch (state) {
      case START:
        if (isalnum(c)) {
          state = ID;
          current_id = c;
        }
        else if (isblank(c)) {
          // ok
        }
        else {
          ERROR("Could not parse reaction side " + str + " (START).");
        }
        break;

      case ID:
        if (isalnum(c)) {
          current_id += c;
        }
        else if (isblank(c) || c == '+') {
          substances.push_back(current_id);
          current_id = "";
          if (isblank(c)) {
            state = AFTER_ID;
          }
          else {
            state = AFTER_PLUS;
          }
        }
        else {
          ERROR("Could not parse reaction side " + str + " (ID).");
        }
        break;

      case AFTER_ID:
        if (c == '+') {
          state = AFTER_PLUS;
        }
        else if (isblank(c)) {
          // ok
        }
        else {
          ERROR("Could not parse reaction side " + str + " (AFTER_ID).");
        }
        break;

      case AFTER_PLUS:
        if (isalnum(c)) {
          state = ID;
          current_id = c;
        }
        else if (isblank(c)) {
          // ok
        }
        else {
          ERROR("Could not parse reaction side " + str + " (AFTER_PLUS).");
        }
        break;
      default:
        assert(false);
    }
  }
  if (current_id != "") {
    substances.push_back(current_id);
  }

  out << "[ ";
  for (size_t i = 0; i < substances.size(); i++) {
    out << substances[i] << "." << NAME_INST << "()";
    print_comma(out, i, substances);
  }
  out << "]";
}

std::string PymcellGenerator::get_filename(const string file_suffix) {
  if (output_files_prefix == "" || output_files_prefix.back() == '/' || output_files_prefix.back() == '\\') {
    return output_files_prefix + file_suffix + PY_EXT;
  }
  else {
    return output_files_prefix + "_" + file_suffix + PY_EXT;
  }
}


std::string PymcellGenerator::get_module_name(const string file_suffix) {
  if (output_files_prefix == "" || output_files_prefix.back() == '/' || output_files_prefix.back() == '\\') {
    return file_suffix;
  }
  else {
    size_t pos = output_files_prefix.find_last_of("/\\");
    if (pos == string::npos) {
      return output_files_prefix + "_" + file_suffix;
    }
    else {
      return output_files_prefix.substr(pos) +  "_" + file_suffix;
    }
  }
}


std::string PymcellGenerator::make_import(const string file_suffix) {
  return "from " + get_module_name(file_suffix) + " import *\n";
}


void PymcellGenerator::open_and_check_file(const string file_suffix, ofstream& out) {
  string file_name = get_filename(file_suffix);
  cout << "Started to generate file " + file_name + ".\n";
  out.open(file_name);
  out.precision(17);
  if (!out.is_open()) {
    throw ConversionError("Could not open file '" + file_name + "' for writing.");
  }
}


void PymcellGenerator::reset() {
  geometry_generated = false;
  instantiation_generated = false;
  observables_generated = false;

  unnamed_rxn_counter = 0;
}


// the aim is to generate as much of output as possible,
// therefore we are using exceptions
bool PymcellGenerator::generate(
    const string& input_file,
    const string& output_files_prefix_
) {
  reset();

  bool failed = false;
  output_files_prefix = output_files_prefix_;

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

  CHECK(generate_parameters(), failed);
  CHECK(generate_subsystem(), failed);
  CHECK(generate_geometry(), failed);
  //CHECK(generate_instantiation(), failed);
  //CHECK(generate_observables(), failed);
  //CHECK(generate_model(), failed);

  return !failed;
}


void PymcellGenerator::generate_parameters() {
  ofstream out;
  open_and_check_file(PARAMETERS_SUFFIX, out);

  out << make_section_comment("simulation setup");

  out << PARAM_SEED << " = 1\n";
  out << PARAM_ITERATIONS << " = " << mcell[KEY_INITIALIZATION][KEY_ITERATIONS].asString() << "\n";
  out << PARAM_TIME_STEP << " = " << mcell[KEY_INITIALIZATION][KEY_TIME_STEP].asString() << "\n";
  out << PARAM_DUMP << " = False\n";

  out.close();
}


void PymcellGenerator::generate_species(ofstream& out) {

  // there must be at least one species
  Value& define_molecules = get_node(mcell, KEY_DEFINE_MOLECULES);
  check_version(KEY_DEFINE_MOLECULES, define_molecules, JSON_DM_VERSION_1638);

  Value& molecule_list = get_node(define_molecules, KEY_MOLECULE_LIST);
  for (Value::ArrayIndex i = 0; i < molecule_list.size(); i++) {
    Value& molecule_list_item = molecule_list[i];
    check_version(KEY_MOLECULE_LIST, molecule_list_item, JSON_DM_VERSION_1632);

    string name = molecule_list_item[KEY_MOL_NAME].asString();
    out << name << " = " << MDOT << NAME_CLASS_SPECIES << "(\n";
    gen_param(out, NAME_NAME, name, true);

    string mol_type = molecule_list_item[KEY_MOL_TYPE].asString();
    CHECK_PROPERTY(mol_type == VALUE_MOL_TYPE_2D || mol_type == VALUE_MOL_TYPE_3D);
    if (mol_type == VALUE_MOL_TYPE_3D) {
      gen_param(out, NAME_DIFFUSION_CONSTANT_3D, molecule_list_item[KEY_DIFFUSION_CONSTANT], false);
    }
    else {
      gen_param(out, NAME_DIFFUSION_CONSTANT_2D, molecule_list_item[KEY_DIFFUSION_CONSTANT], false);
    }
    out << CTOR_END;
  }
}


void PymcellGenerator::generate_reaction_rules(ofstream& out) {

  if (!mcell.isMember(KEY_DEFINE_REACTIONS)) {
    return;
  }

  Value& define_reactions = get_node(mcell, KEY_DEFINE_REACTIONS);
  check_version(KEY_DEFINE_REACTIONS, define_reactions, JSON_DM_VERSION_1638);

  Value& reaction_list = get_node(define_reactions, KEY_REACTION_LIST);
  for (Value::ArrayIndex i = 0; i < reaction_list.size(); i++) {
    Value& reaction_list_item = reaction_list[i];
    check_version(KEY_MOLECULE_LIST, reaction_list_item, JSON_DM_VERSION_1330);

    string name = reaction_list_item[KEY_RXN_NAME].asString();
    if (name == "") {
      name = UNNAMED_REACTION_RULE_PREFIX + to_string(unnamed_rxn_counter);
      unnamed_rxn_counter++;
    }
    out << name << " = " << MDOT << NAME_CLASS_SPECIES << "(\n";
    gen_param(out, NAME_NAME, name, true);

    // single line for now
    out << IND << NAME_REACTANTS << " = ";
    gen_rxn_substance_inst(out, reaction_list_item[KEY_REACTANTS]);
    out << ",\n";

    out << IND << NAME_PRODUCTS << " = ";
    gen_rxn_substance_inst(out, reaction_list_item[KEY_PRODUCTS]);
    out << ",\n";

    CHECK_PROPERTY(reaction_list_item[KEY_VARIABLE_RATE_SWITCH].asBool() == false && "Not supported yet");

    string rxn_type = reaction_list_item[KEY_RXN_TYPE].asString();
    CHECK_PROPERTY(rxn_type == VALUE_IRREVERSIBLE || rxn_type == VALUE_REVERSIBLE);
    bool is_reversible = rxn_type == VALUE_REVERSIBLE;

    gen_param(out, NAME_FWD_RATE, reaction_list_item[KEY_FWD_RATE], is_reversible);
    if (is_reversible) {
      gen_param(out, NAME_FWD_RATE, reaction_list_item[KEY_BKWD_RATE], false);
    }

    out << CTOR_END;
  }

}


void PymcellGenerator::generate_subsystem() {
  ofstream out;
  open_and_check_file(SUBSYSTEM_SUFFIX, out);

  out << MCELL_IMPORT;
  out << make_import(PARAMETERS_SUFFIX);
  out << "\n";
  out << make_section_comment(SUBSYSTEM_SUFFIX);

  generate_species(out);
  generate_reaction_rules(out);

  out.close();
}


void PymcellGenerator::generate_single_geometry_object(
    ostream& out, const int index, Value& object) {

  string parent_name = S(KEY_OBJECT_LIST) + "[" + to_string(index) + "]";

  string name = get_node(parent_name, object, KEY_NAME).asString();
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
      print_comma(out, k, vertex);
    }
    out << "]";
    print_comma(out, i, vertex_list);
    out << "\n";
  }
  out << "] # " << id_vertex_list << "\n\n";

  // element_connections
  string id_element_connections = name + "_" + NAME_ELEMENT_CONNECTIONS;
  out << id_element_connections << " = [\n";
  for (Value::ArrayIndex i = 0; i < element_connections.size(); i++) {
    out << IND << "[";
    Value& element = element_connections[i];
    for (Value::ArrayIndex k = 0; k < element.size(); k++) {
      out << element[k].asInt();
      print_comma(out, k, element);
    }
    out << "]";
    print_comma(out, i, element_connections);
    out << "\n";
  }
  out << "] # " << id_element_connections << "\n\n";

  // surface regions
  vector<string> sr_global_names;
  if (object.isMember(KEY_DEFINE_SURFACE_REGIONS)) {
    Value& define_surface_regions = object[KEY_DEFINE_SURFACE_REGIONS];
    for (Value::ArrayIndex i = 0; i < define_surface_regions.size(); i++) {

      string sr_name = get_node(
          KEY_DEFINE_SURFACE_REGIONS, define_surface_regions[i], KEY_NAME).asString();
      Value& include_elements = get_node(
          KEY_DEFINE_SURFACE_REGIONS, define_surface_regions[i], KEY_INCLUDE_ELEMENTS);

      string sr_global_name = name + "_" + sr_name;
      string sr_element_connections_name = sr_global_name + "_" + NAME_ELEMENT_CONNECTIONS;
      out << sr_element_connections_name << " = [";
      for (Value::ArrayIndex k = 0; k < include_elements.size(); k++) {

        if (k % 16 == 0) {
          out << "\n" << IND;
        }
        out << include_elements[k].asInt();
        print_comma(out, k, include_elements);
      }
      out << "\n] #" << sr_element_connections_name << "\n\n";

      out << sr_global_name << " = " << MDOT << NAME_CLASS_SURFACE_REGION << "(\n";
      out << IND << NAME_NAME << " = '" << name << "',\n";
      out << IND << NAME_ELEMENT_CONNECTIONS << " = " << sr_element_connections_name << "\n";
      out << ")\n\n";

      sr_global_names.push_back(sr_global_name);
    }
  }

  // object creation itself
  out << name << " = " << MDOT << NAME_CLASS_GEOMETRY_OBJECT << "(\n";
  out << IND << NAME_NAME << " = '" << name << "',\n";
  out << IND << NAME_VERTEX_LIST << " = " << id_vertex_list << ",\n";
  out << IND << NAME_ELEMENT_CONNECTIONS << " = " << id_element_connections << ",\n";
  out << IND << NAME_SURFACE_REGIONS << " = [";
  for (size_t i = 0; i < sr_global_names.size(); i++) {
    out << sr_global_names[i];
    print_comma(out, i, sr_global_names);
  }
  out << "]\n)\n";

  out << make_end_block_comment(name);
}


void PymcellGenerator::generate_geometry() {

  // TODO: check versions

  if (!mcell.isMember(KEY_GEOMETRICAL_OBJECTS)) {
    return;
  }
  Value& geometrical_objects = get_node(mcell, KEY_GEOMETRICAL_OBJECTS);
  if (!geometrical_objects.isMember(KEY_OBJECT_LIST)) {
    return;
  }

  ofstream out;
  open_and_check_file(GEOMETRY_SUFFIX, out);

  out << MCELL_IMPORT;

  Value& object_list = get_node(geometrical_objects, KEY_OBJECT_LIST);
  for (Value::ArrayIndex i = 0; i < object_list.size(); i++) {
    Value& object = object_list[i];
    generate_single_geometry_object(out, i, object);
  }

  out.close();
  geometry_generated = true;
}


} /* namespace MCell */
