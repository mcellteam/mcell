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

#include <regex>

#include "pymcell_generator.h"

#include "generator_utils.h"

using namespace std;
using namespace MCell::API;

namespace MCell {

using Json::Value;


string PymcellGenerator::get_filename(const string file_suffix) {
  if (output_files_prefix == "" || output_files_prefix.back() == '/' || output_files_prefix.back() == '\\') {
    return output_files_prefix + file_suffix + PY_EXT;
  }
  else {
    return output_files_prefix + "_" + file_suffix + PY_EXT;
  }
}


string PymcellGenerator::get_module_name(const string file_suffix) {
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


string PymcellGenerator::make_import(const string file_suffix) {
  return "from " + get_module_name(file_suffix) + " import *\n";
}


void PymcellGenerator::open_and_check_file(const string file_suffix, ofstream& out) {
  string file_name = get_filename(file_suffix);
  cout << "Generating file " + file_name + ".\n";
  out.open(file_name);
  out.precision(17);
  if (!out.is_open()) {
    throw ConversionError("Could not open file '" + file_name + "' for writing.");
  }
}


void PymcellGenerator::reset() {
  geometry_generated = false;
  observables_generated = false;
  unnamed_rxn_counter = 0;
  all_species_names.clear();
}


// the aim is to generate as much of output as possible,
// therefore we are using exceptions
bool PymcellGenerator::generate(
    const string& input_file,
    const string& output_files_prefix_,
    const bool debug_mode_
) {
  reset();

  bool failed = false;
  output_files_prefix = output_files_prefix_;
  debug_mode = true;

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
  std::vector<std::string> geometry_names;
  CHECK(geometry_names = generate_geometry(), failed);
  CHECK(generate_instantiation(geometry_names), failed);
  CHECK(generate_observables(), failed);
  CHECK(generate_model(), failed);

  return !failed;
}


void PymcellGenerator::generate_parameters() {
  ofstream out;
  open_and_check_file(PARAMETERS, out);

  out << make_section_comment("simulation setup");

  out << PARAM_SEED << " = 1\n";
  out << PARAM_ITERATIONS << " = " << mcell[KEY_INITIALIZATION][KEY_ITERATIONS].asString() << "\n";
  out << PARAM_TIME_STEP << " = " << mcell[KEY_INITIALIZATION][KEY_TIME_STEP].asString() << "\n";
  out << PARAM_DUMP << " = " << (debug_mode ? "True" : "False") << "\n";

  out.close();
}


vector<string> PymcellGenerator::generate_species(ofstream& out) {
  vector<string> species_names;

  // there must be at least one species
  Value& define_molecules = get_node(mcell, KEY_DEFINE_MOLECULES);
  check_version(KEY_DEFINE_MOLECULES, define_molecules, JSON_DM_VERSION_1638);

  Value& molecule_list = get_node(define_molecules, KEY_MOLECULE_LIST);
  for (Value::ArrayIndex i = 0; i < molecule_list.size(); i++) {
    Value& molecule_list_item = molecule_list[i];
    check_version(KEY_MOLECULE_LIST, molecule_list_item, JSON_DM_VERSION_1632);

    string name = molecule_list_item[KEY_MOL_NAME].asString();
    species_names.push_back(name);
    gen_ctor_call(out, name, NAME_CLASS_SPECIES);
    gen_param(out, NAME_NAME, name, true);

    string mol_type = molecule_list_item[KEY_MOL_TYPE].asString();
    CHECK_PROPERTY(mol_type == VALUE_MOL_TYPE_2D || mol_type == VALUE_MOL_TYPE_3D);
    if (mol_type == VALUE_MOL_TYPE_3D) {
      gen_param_double(out, NAME_DIFFUSION_CONSTANT_3D, molecule_list_item[KEY_DIFFUSION_CONSTANT], false);
    }
    else {
      gen_param_double(out, NAME_DIFFUSION_CONSTANT_2D, molecule_list_item[KEY_DIFFUSION_CONSTANT], false);
    }
    out << CTOR_END;
  }

  return species_names;
}


vector<string> PymcellGenerator::generate_reaction_rules(ofstream& out) {
  vector<string> rxn_names;

  if (!mcell.isMember(KEY_DEFINE_REACTIONS)) {
    return rxn_names;
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
    rxn_names.push_back(name);
    gen_ctor_call(out, name, NAME_CLASS_REACTION_RULE);
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

    gen_param_double(out, NAME_FWD_RATE, reaction_list_item[KEY_FWD_RATE], is_reversible);
    if (is_reversible) {
      gen_param_double(out, NAME_FWD_RATE, reaction_list_item[KEY_BKWD_RATE], false);
    }

    out << CTOR_END;
  }

  return rxn_names;
}


void PymcellGenerator::generate_subsystem() {
  ofstream out;
  open_and_check_file(SUBSYSTEM, out);

  out << MCELL_IMPORT;
  out << make_import(PARAMETERS);
  out << "\n";
  out << make_section_comment(SUBSYSTEM);

  all_species_names = generate_species(out);
  vector<string> rxn_names = generate_reaction_rules(out);

  gen_ctor_call(out, SUBSYSTEM, NAME_CLASS_SUBSYSTEM, false);
  for (string& s: all_species_names) {
    gen_method_call(out, SUBSYSTEM, NAME_ADD_SPECIES, s);
  }
  for (string& r: rxn_names) {
    gen_method_call(out, SUBSYSTEM, NAME_ADD_REACTION_RULE, r);
  }

  out.close();
}


string PymcellGenerator::generate_single_geometry_object(
    ofstream& out, const int index, Value& object) {

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
      gen_comma(out, k, vertex);
    }
    out << "]";
    gen_comma(out, i, vertex_list);
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

      string sr_name = get_node(
          KEY_DEFINE_SURFACE_REGIONS, define_surface_regions[i], KEY_NAME).asString();
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
      out << ")\n\n";

      sr_global_names.push_back(sr_global_name);
    }
  }

  // object creation itself
  out << name << " = " << MDOT << NAME_CLASS_GEOMETRY_OBJECT << "(\n";
  gen_param(out, NAME_NAME, name, true);
  gen_param_id(out, NAME_VERTEX_LIST, id_vertex_list, true);
  gen_param_id(out, NAME_ELEMENT_CONNECTIONS, id_element_connections, true);
  out << IND << NAME_SURFACE_REGIONS << " = [";
  for (size_t i = 0; i < sr_global_names.size(); i++) {
    out << sr_global_names[i];
    print_comma(out, i, sr_global_names);
  }
  out << "]\n)\n";

  out << make_end_block_comment(name);

  return name;
}


vector<string> PymcellGenerator::generate_geometry() {

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

  Value& object_list = get_node(geometrical_objects, KEY_OBJECT_LIST);
  for (Value::ArrayIndex i = 0; i < object_list.size(); i++) {
    Value& object = object_list[i];
    string name = generate_single_geometry_object(out, i, object);
    geometry_objects.push_back(name);
  }

  out.close();
  geometry_generated = true;

  return geometry_objects;
}


void gen_region_expr_assignment_for_rel_site(ofstream& out, string region_expr) {
  // we don't really need to parse the expression, just dump it in a right way
  // example: Cell[ALL] - (Organelle_1[ALL] + Organelle_2[ALL])
  // remove [ALL]
  // replace [text] with _text

  regex pattern_all("\\[ALL\\]");
  region_expr = regex_replace(region_expr, pattern_all, "");

  regex pattern_surf_region("\\[(^\\])\\]");
  region_expr = regex_replace(region_expr, pattern_surf_region, "_$1");

  gen_param_id(out, NAME_REGION, region_expr, true);
}


vector<string> PymcellGenerator::generate_release_sites(ofstream& out) {
  vector<string> release_site_names;

  if (!mcell.isMember(KEY_RELEASE_SITES)) {
    return release_site_names;
  }

  Json::Value& release_sites = mcell[KEY_RELEASE_SITES];
  check_version(KEY_RELEASE_SITES, release_sites, JSON_DM_VERSION_1638);
  Json::Value& release_site_list = release_sites[KEY_RELEASE_SITE_LIST];
  for (Value::ArrayIndex i = 0; i < release_site_list.size(); i++) {
    Value& release_site_item = release_site_list[i];
    check_version(KEY_RELEASE_SITE_LIST, release_site_item, JSON_DM_VERSION_1330);
    string name = release_site_item[KEY_NAME].asString();
    release_site_names.push_back(name);

    gen_ctor_call(out, name, NAME_CLASS_RELEASE_SITE);
    gen_param(out, NAME_NAME, name, true);
    gen_param_id(out, NAME_SPECIES, release_site_item[KEY_MOLECULE], true);

    string shape = release_site_item[KEY_SHAPE].asString();
    if (shape == VALUE_SPHERICAL) {
      gen_param_enum(out, NAME_SHAPE, NAME_ENUM_SHAPE, NAME_EV_SPHERICAL, true);
      gen_param_vec3(
          out, NAME_LOCATION,
          release_site_item[KEY_LOCATION_X], release_site_item[KEY_LOCATION_Y], release_site_item[KEY_LOCATION_Z],
          true
      );
    }
    else if (shape == VALUE_OBJECT) {
      gen_region_expr_assignment_for_rel_site(out, release_site_item[KEY_OBJECT_EXPR].asString());
    }
    else {
      ERROR("Shape " + shape + " is not supported yet");
    }

    string quantity_type = release_site_item[KEY_QUANTITY_TYPE].asString();
    if (quantity_type == VALUE_NUMBER_TO_RELEASE) {
      gen_param_int(out, NAME_NUMBER_TO_RELEASE, release_site_item[KEY_QUANTITY], false);
    }
    else {
      ERROR("Quantity type " + quantity_type + " is not supported yet");
    }

    out << ")\n\n";
  }

  return release_site_names;
}


void PymcellGenerator::generate_instantiation(const vector<string>& geometry_objects) {

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

  vector<string> release_sites = generate_release_sites(out);

  gen_ctor_call(out, INSTANTIATION, NAME_CLASS_INSTANTIATION_DATA, false);
  for (const string& s: geometry_objects) {
    gen_method_call(out, INSTANTIATION, NAME_ADD_GEOMETRY_OBJECT, s);
  }
  for (const string& r: release_sites) {
    gen_method_call(out, INSTANTIATION, NAME_ADD_RELEASE_SITE, r);
  }

  out.close();
}


vector<string> PymcellGenerator::get_species_to_visualize() {
  vector<string> res;

  Value& define_molecules = get_node(mcell, KEY_DEFINE_MOLECULES);
  check_version(KEY_DEFINE_MOLECULES, define_molecules, JSON_DM_VERSION_1638);

  Value& molecule_list = get_node(define_molecules, KEY_MOLECULE_LIST);
  for (Value::ArrayIndex i = 0; i < molecule_list.size(); i++) {
    Value& molecule_list_item = molecule_list[i];
    check_version(KEY_MOLECULE_LIST, molecule_list_item, JSON_DM_VERSION_1632);

    if (molecule_list_item[KEY_EXPORT_VIZ].asBool()) {
      res.push_back(molecule_list_item[KEY_MOL_NAME].asString());
    }
  }

  return res;
}


vector<string> PymcellGenerator::generate_viz_outputs(ofstream& out) {
  vector<string> viz_output_names;

  if (!mcell.isMember(KEY_VIZ_OUTPUT)) {
    return viz_output_names;
  }

  Json::Value& viz_output = mcell[KEY_VIZ_OUTPUT];
  check_version(KEY_VIZ_OUTPUT, viz_output, JSON_DM_VERSION_1638);

  string name = VIZ_OUTPUT_NAME; // there is only one in datamodel now
  viz_output_names.push_back(name);

  // CHECK_PROPERTY(viz_output[KEY_ALL_ITERATIONS].asBool()); // don't care
  CHECK_PROPERTY(viz_output[KEY_START].asString() == "0");

  gen_ctor_call(out, name, NAME_CLASS_VIZ_OUTPUT);

  // mode is ascii by default, this information is not in datamodel (AFAIK)
  gen_param_enum(out, NAME_MODE, NAME_ENUM_VIZ_MODE, NAME_EV_ASCII, true);
  gen_param(out, NAME_FILENAME_PREFIX, DEFAULT_VIZ_OUTPUT_FILENAME_PREFIX, true);

  // species_list
  std::vector<std::string> viz_species;
  if (viz_output[KEY_EXPORT_ALL].asBool()) {
    viz_species = all_species_names;
  }
  else {
    viz_species = get_species_to_visualize();
  }
  gen_param_list(out, NAME_SPECIES_LIST, viz_species, true);

  gen_param_int(out, NAME_EVERY_N_TIMESTEPS, viz_output[KEY_STEP], false);

  // ignoring KEY_END
  out << ")\n\n";

  return viz_output_names;
}


vector<string> PymcellGenerator::generate_counts(ofstream& out) {
  vector<string> counts;
  // TODO

  return counts;
}


void PymcellGenerator::generate_observables() {

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

  vector<string> viz_outputs = generate_viz_outputs(out);

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


void PymcellGenerator::generate_config(ofstream& out) {
  out << make_section_comment("configuration");

  // using values from generated parameters.py
  gen_assign(out, MODEL, NAME_CONFIG, NAME_TIME_STEP, PARAM_TIME_STEP);
  gen_assign(out, MODEL, NAME_CONFIG, NAME_SEED, PARAM_SEED);
  gen_assign(out, MODEL, NAME_CONFIG, NAME_TOTAL_ITERATIONS_HINT, PARAM_ITERATIONS);
  out << "\n";

  Json::Value& initialization = mcell[KEY_INITIALIZATION];
  Json::Value& partitions = initialization[KEY_PARTITIONS];

  // check that all the x,y,z values are the same (comparing strings)
  CHECK_PROPERTY(partitions[KEY_X_START].asString() == partitions[KEY_Y_START].asString());
  CHECK_PROPERTY(partitions[KEY_Y_START].asString() == partitions[KEY_Z_START].asString());
  CHECK_PROPERTY(partitions[KEY_X_END].asString() == partitions[KEY_Y_END].asString());
  CHECK_PROPERTY(partitions[KEY_Y_END].asString() == partitions[KEY_Z_END].asString());
  CHECK_PROPERTY(partitions[KEY_X_STEP].asString() == partitions[KEY_Y_STEP].asString());
  CHECK_PROPERTY(partitions[KEY_Z_STEP].asString() == partitions[KEY_Z_STEP].asString());

  double x_start = stod(partitions[KEY_X_START].asString());
  double x_end = stod(partitions[KEY_X_END].asString());
  double x_step = stod(partitions[KEY_X_STEP].asString());
  CHECK_PROPERTY(fabs(x_start + x_end) < 1e-12 && "-x_start must be equal to x_end");

  gen_assign(out, MODEL, NAME_CONFIG, NAME_PARTITION_DIMENSION, x_end * 2);
  gen_assign(out, MODEL, NAME_CONFIG, NAME_SUBPARTITION_DIMENSION, x_step);
}


void PymcellGenerator::generate_model() {
  ofstream out;
  open_and_check_file(MODEL, out);

  out << INTERPRETER;
  out << BASE_MODEL_IMPORTS;
  out << MCELL_DIR_SETUP;

  out << MCELL_IMPORT;

  out << make_import(PARAMETERS);
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

  gen_method_call(out, MODEL, NAME_RUN_ITERATIONS, PARAM_ITERATIONS);
  gen_method_call(out, MODEL, NAME_END_SIMULATION);
}


} // namespace MCell
