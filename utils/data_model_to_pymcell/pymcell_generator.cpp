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
#include <algorithm>

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
  out.precision(FLOAT_OUT_PRECISION);
  if (!out.is_open()) {
    throw ConversionError("Could not open file '" + file_name + "' for writing.");
  }
}


void PymcellGenerator::reset() {
  geometry_generated = false;
  observables_generated = false;
  unnamed_rxn_counter = 0;
  unnamed_surf_class_counter = 0;
  all_species_names.clear();
  all_reaction_rules_names.clear();
  all_count_term_names.clear();
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


void PymcellGenerator::get_surface_class_property_info(
    Value& property,
    string& name, string& type_name, string& affected_mols, string& orientation
) {
  check_version(KEY_SURFACE_CLASS_PROP_LIST, property, JSON_DM_VERSION_1756);
  CHECK_PROPERTY(property[KEY_CLAMP_VALUE] = "0");

  affected_mols = property[KEY_AFFECTED_MOLS].asString();
  if (affected_mols == ALL_MOLECULES) {
    affected_mols = S(MDOT) + NAME_CV_AllMolecules;
  }
  else if (affected_mols == ALL_VOLUME_MOLECULES) {
    affected_mols = S(MDOT) + NAME_CV_AllVolumeMolecules;
  }
  else if (affected_mols == ALL_SURFACE_MOLECULES) {
    affected_mols = S(MDOT) + NAME_CV_AllSurfaceMolecules;
  }
  else if (affected_mols == VALUE_SINGLE) {
    affected_mols = property[KEY_MOLECULE].asString();
  }

  orientation = convert_orientation(property[KEY_NAME].asString());

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
  else {
    ERROR(S("Invalid ") + KEY_SURF_CLASS_TYPE + " " + surf_class_type + ".");
  }

  name = property[KEY_NAME].asString();
  if (name == "") {
    string tn_lower = type_name;
    transform(tn_lower.begin(), tn_lower.end(), tn_lower.begin(), ::tolower);
    // to be sure that there is no conflict, let's put a number at the end
    name = SURFACE_CLASS_PREFIX + tn_lower + "_" + affected_mols + to_string(unnamed_surf_class_counter);
    unnamed_surf_class_counter++;
  }
}


vector<string> PymcellGenerator::generate_surface_classes(ofstream& out) {
  vector<string> sc_names;

  if (!mcell.isMember(KEY_DEFINE_SURFACE_CLASSES)) {
    return sc_names;
  }

  Value& define_surface_classes = get_node(mcell, KEY_DEFINE_SURFACE_CLASSES);
  check_version(KEY_DEFINE_SURFACE_CLASSES, define_surface_classes, JSON_DM_VERSION_1638);
  Value& surface_class_list = get_node(define_surface_classes, KEY_SURFACE_CLASS_LIST);

  for (Value::ArrayIndex i = 0; i < surface_class_list.size(); i++) {
    Value& surface_class_list_item = surface_class_list[i];

    Value& surface_class_prop_list = get_node(surface_class_list_item, KEY_SURFACE_CLASS_PROP_LIST);

    vector<string> sc_prop_names;
    if (surface_class_prop_list.size() > 1) {
      for (Value::ArrayIndex i = 0; i < surface_class_prop_list.size(); i++) {
        Value& surface_class_prop_item = surface_class_prop_list[i];

        string name, type_name, affected_mols, orientation_name;
        get_surface_class_property_info(surface_class_prop_item, name, type_name, affected_mols, orientation_name);

        sc_prop_names.push_back(name);
        gen_ctor_call(out, name, NAME_CLASS_SURFACE_PROPERTY, true);
        gen_param_enum(out, NAME_TYPE, NAME_ENUM_SURFACE_PROPERTY_TYPE, type_name, true);
        gen_param_id(out, NAME_AFFECTED_SPECIES, affected_mols, orientation_name != "");
        if (orientation_name != "") {
          gen_param_id(out, NAME_ORIENTATION, orientation_name, false);
        }
        out << CTOR_END;
      }
    }

    string name = surface_class_list_item[KEY_NAME].asString();
    sc_names.push_back(name);
    gen_ctor_call(out, name, NAME_CLASS_SURFACE_CLASS, true);
    gen_param(out, NAME_NAME, name, true);

    if (!sc_prop_names.empty()) {
      // use a list of properties
      gen_param_list(out, NAME_PROPERTIES, sc_prop_names, false);
    }
    else {
      // simplified setup, directly set members
      string name, type_name, affected_mols, orientation_name;
      get_surface_class_property_info(surface_class_prop_list[0], name, type_name, affected_mols, orientation_name);

      gen_param_enum(out, NAME_TYPE, NAME_ENUM_SURFACE_PROPERTY_TYPE, type_name, true);
      gen_param_id(out, NAME_AFFECTED_SPECIES, affected_mols, orientation_name != "");
      if (orientation_name != "") {
        gen_param_id(out, NAME_ORIENTATION, orientation_name, false);
      }
    }
    out << CTOR_END;
  }

  return sc_names;
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
  vector<string> surface_class_names = generate_surface_classes(out);
  all_reaction_rules_names = generate_reaction_rules(out);

  gen_ctor_call(out, SUBSYSTEM, NAME_CLASS_SUBSYSTEM, false);
  for (string& s: all_species_names) {
    gen_method_call(out, SUBSYSTEM, NAME_ADD_SPECIES, s);
  }
  for (string& sc: surface_class_names) {
    gen_method_call(out, SUBSYSTEM, NAME_ADD_SURFACE_CLASS, sc);
  }
  for (string& r: all_reaction_rules_names) {
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

  regex pattern_surf_region("\\[([^\\]]*)\\]");
  region_expr = regex_replace(region_expr, pattern_surf_region, "_$1");

  gen_param_id(out, NAME_REGION, region_expr, true);
}


vector<string> PymcellGenerator::generate_release_sites(ofstream& out) {
  vector<string> release_site_names;

  if (!mcell.isMember(KEY_RELEASE_SITES)) {
    return release_site_names;
  }

  Value& release_sites = mcell[KEY_RELEASE_SITES];
  check_version(KEY_RELEASE_SITES, release_sites, JSON_DM_VERSION_1638);
  Value& release_site_list = release_sites[KEY_RELEASE_SITE_LIST];
  for (Value::ArrayIndex i = 0; i < release_site_list.size(); i++) {
    Value& release_site_item = release_site_list[i];
    check_version(KEY_RELEASE_SITE_LIST, release_site_item, JSON_DM_VERSION_1330);
    string name = release_site_item[KEY_NAME].asString();
    release_site_names.push_back(name);

    gen_ctor_call(out, name, NAME_CLASS_RELEASE_SITE);
    gen_param(out, NAME_NAME, name, true);
    gen_param_id(out, NAME_SPECIES, release_site_item[KEY_MOLECULE], true);

    string orientation = convert_orientation(release_site_item[KEY_ORIENT].asString());
    if (orientation != "") {
      gen_param_enum(out, NAME_ORIENTATION, NAME_ENUM_ORIENTATION, orientation, true);
    }

    string shape = release_site_item[KEY_SHAPE].asString();
    if (shape == VALUE_SPHERICAL) {
      gen_param_enum(out, NAME_SHAPE, NAME_ENUM_SHAPE, NAME_EV_SPHERICAL, true);
      gen_param_vec3(
          out, NAME_LOCATION,
          release_site_item[KEY_LOCATION_X], release_site_item[KEY_LOCATION_Y], release_site_item[KEY_LOCATION_Z],
          true
      );
      gen_param_double(out, NAME_SITE_DIAMETER, release_site_item[KEY_SITE_DIAMETER], true);
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


void PymcellGenerator::generate_surface_classes_assignment(ofstream& out) {
  if (!mcell.isMember(KEY_MODIFY_SURFACE_REGIONS)) {
    return;
  }

  Value& modify_surface_regions = mcell[KEY_MODIFY_SURFACE_REGIONS];
  check_version(KEY_MODIFY_SURFACE_REGIONS, modify_surface_regions, JSON_DM_VERSION_1638);
  Value& modify_surface_regions_list = modify_surface_regions[KEY_MODIFY_SURFACE_REGIONS_LIST];
  for (Value::ArrayIndex i = 0; i < modify_surface_regions_list.size(); i++) {
    Value& modify_surface_regions_item = modify_surface_regions_list[i];
    check_version(KEY_MODIFY_SURFACE_REGIONS_LIST, modify_surface_regions_item, JSON_DM_VERSION_1330);

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

    string surf_class_name = modify_surface_regions_item[KEY_SURF_CLASS_NAME].asString();

    // generate the assignment
    out << obj_or_region_name << "." << NAME_SURFACE_CLASS << " = " << surf_class_name << "\n";
  }
  out << "\n";
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
  out << make_section_comment("release sites");

  vector<string> release_sites = generate_release_sites(out);

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

  Value& viz_output = mcell[KEY_VIZ_OUTPUT];
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


void PymcellGenerator::process_single_count_term(
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

  if (find(all_species_names.begin(), all_species_names.end(), what_to_count) != all_species_names.end()) {
    rxn_not_mol = false;
  }
  else if (find(all_reaction_rules_names.begin(), all_reaction_rules_names.end(), what_to_count) != all_reaction_rules_names.end()) {
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
}


string PymcellGenerator::generate_count_terms_for_expression(ofstream& out, const string& mdl_string) {
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

      out << ")\n\n";
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


vector<string> PymcellGenerator::generate_counts(ofstream& out) {
  vector<string> counts;

  if (!mcell.isMember(KEY_REACTION_DATA_OUTPUT)) {
    return counts;
  }

  Value& reaction_data_output = get_node(mcell, KEY_REACTION_DATA_OUTPUT);
  check_version(KEY_DEFINE_MOLECULES, reaction_data_output, JSON_DM_VERSION_1800);

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

    string rxn_or_mol = reaction_output_item[KEY_RXN_OR_MOL].asString();
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
    }
    else if (rxn_or_mol == VALUE_REACTION) {
      single_term = true;
      rxn_not_mol = true;
      what_to_count = reaction_output_item[KEY_REACTION_NAME].asString();
      where_to_count = reaction_output_item[KEY_REGION_NAME].asString();
    }
    else if (rxn_or_mol == VALUE_MOLECULE) {
      single_term = true;
      rxn_not_mol = false;
      what_to_count = reaction_output_item[KEY_MOLECULE_NAME].asString();
      where_to_count = reaction_output_item[KEY_REGION_NAME].asString();
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

    gen_param_id(out, NAME_FILENAME,
        DEFAULT_RXN_OUTPUT_FILENAME_PREFIX + mdl_file_prefix + ".dat'", rxn_step != "");

    if (rxn_step != "") {
      gen_param_int(out, NAME_EVERY_N_TIMESTEPS, rxn_step, false);
    }

    out << ")\n\n";
  }

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

  Value& initialization = mcell[KEY_INITIALIZATION];
  Value& partitions = initialization[KEY_PARTITIONS];

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
