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

#include <fstream>

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


void PythonGenerator::get_surface_class_property_info(
    Value& property,
    string& name, string& type_name, string& affected_mols, string& orientation
) {
  check_version(KEY_SURFACE_CLASS_PROP_LIST, property, VER_DM_2015_11_08_1756);
  CHECK_PROPERTY(property[KEY_CLAMP_VALUE] = "0");

  affected_mols = property[KEY_AFFECTED_MOLS].asString();
  if (affected_mols == ::ALL_MOLECULES) {
    affected_mols = S(MDOT) + NAME_CV_AllMolecules;
  }
  else if (affected_mols == ::ALL_VOLUME_MOLECULES) {
    affected_mols = S(MDOT) + NAME_CV_AllVolumeMolecules;
  }
  else if (affected_mols == ::ALL_SURFACE_MOLECULES) {
    affected_mols = S(MDOT) + NAME_CV_AllSurfaceMolecules;
  }
  else if (affected_mols == VALUE_SINGLE) {
    affected_mols = property[KEY_MOLECULE].asString();
  }

  orientation = convert_orientation(property[KEY_SURF_CLASS_ORIENT].asString());

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

    string name = make_id(surface_class_list_item[KEY_NAME].asString());
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
}


void PythonGenerator::generate_variable_rate(const std::string& rate_array_name, Json::Value& variable_rate_text) {
  // append to the parameters file
  ofstream out;
  open_and_check_file_w_prefix(output_files_prefix, PARAMETERS, out, true);

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


std::string PythonGenerator::generate_single_reaction_rule(std::ostream& out, Json::Value& reaction_list_item) {
  check_version(KEY_MOLECULE_LIST, reaction_list_item, VER_DM_2018_01_11_1330);

  // TODO: BNG rules support

  string name = reaction_list_item[KEY_RXN_NAME].asString();
  if (name == "") {
    bool ok = convert_reaction_name(reaction_list_item[KEY_NAME].asString(), name);

    if (!ok) {
      name = UNNAMED_REACTION_RULE_PREFIX + to_string(unnamed_rxn_counter);
      unnamed_rxn_counter++;
    }
  }
  gen_ctor_call(out, name, NAME_CLASS_REACTION_RULE);
  gen_param(out, NAME_NAME, name, true);

  // single line for now
  out << IND << NAME_REACTANTS << " = ";
  gen_rxn_substance_inst(out, reaction_list_item[KEY_REACTANTS]);
  out << ",\n";

  out << IND << NAME_PRODUCTS << " = ";
  gen_rxn_substance_inst(out, reaction_list_item[KEY_PRODUCTS]);
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


void PythonGenerator::generate_geometry(std::ostream& out, std::vector<std::string>& geometry_objects) {

  Value& geometrical_objects = get_node(mcell, KEY_GEOMETRICAL_OBJECTS);
  if (!geometrical_objects.isMember(KEY_OBJECT_LIST)) {
    return;
  }
  Value& object_list = get_node(geometrical_objects, KEY_OBJECT_LIST);
  for (Value::ArrayIndex i = 0; i < object_list.size(); i++) {
    Value& object = object_list[i];
    string name = generate_single_geometry_object(out, i, object);
    geometry_objects.push_back(name);
  }
}


bool PythonGenerator::is_volume_mol_type(const std::string& mol_type_name) {

  Value& define_molecules = get_node(mcell, KEY_DEFINE_MOLECULES);
  check_version(KEY_DEFINE_MOLECULES, define_molecules, VER_DM_2014_10_24_1638);

  Value& molecule_list = get_node(define_molecules, KEY_MOLECULE_LIST);
  for (Value::ArrayIndex i = 0; i < molecule_list.size(); i++) {
    Value& molecule_list_item = molecule_list[i];
    check_version(KEY_MOLECULE_LIST, molecule_list_item, VER_DM_2018_10_16_1632);

    string name = molecule_list_item[KEY_MOL_NAME].asString();
    if (name != mol_type_name) {
      continue;
    }

    string mol_type = molecule_list_item[KEY_MOL_TYPE].asString();
    CHECK_PROPERTY(mol_type == VALUE_MOL_TYPE_2D || mol_type == VALUE_MOL_TYPE_3D);
    return mol_type == VALUE_MOL_TYPE_3D;
  }

  ERROR("Could not find species or molecule type " + mol_type_name + ".");
}


static bool is_simple_species(const std::string& species_name) {
  return species_name.find('(') == string::npos;
}


static void get_mol_types_in_species(const std::string& species_name, vector<string>& mol_types) {
  mol_types.clear();

  size_t i = 0;
  string current_name;
  bool in_name = true;
  while (i < species_name.size()) {
    char c = species_name[i];
    if (c == '(') {
      in_name = false;
      assert(current_name != "");
      mol_types.push_back(current_name);
      current_name = "";
    }
    else if (c == '.') {
      in_name = true;
    }
    else if (in_name && !isspace(c)) {
      current_name += c;
    }
    i++;
  }

  if (current_name != "") {
    mol_types.push_back(current_name);
  }
}


bool PythonGenerator::is_volume_species(const std::string& species_name) {
  if (is_simple_species(species_name)) {
    return is_volume_mol_type(species_name);
  }
  else {
    vector<string> mol_types;
    get_mol_types_in_species(species_name, mol_types);
    for (string& mt: mol_types) {
      if (!is_volume_mol_type(mt)) {
        // surface
        return false;
      }
    }
    return true;
  }
}


string remove_compartments(const std::string& species_name) {
  size_t i = 0;
  string res;
  bool in_compartment = false;
  while (i < species_name.size()) {
    char c = species_name[i];
    if (c == '@') {
      assert(!in_compartment);
      in_compartment = true;
    }
    else if (in_compartment && (!isalnum(c) && c != '_')) {
      in_compartment = false;
    }
    else if (!in_compartment) {
      res += c;
    }
    i++;
  }
  return res;
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
      out << "    " << MDOT << NAME_CLASS_MOLECULE_RELEASE_INFO << "(";

      string species_name = release_site_item[KEY_MOLECULE].asString();
      if (!bng_mode) {
        gen_param_id(out, NAME_SPECIES, species_name, true);
      }
      else {
        gen_param(out, NAME_BNGL_SPECIES, remove_compartments(species_name), true);
      }

      Value& point = points_list[i];
      if (point.size() != 3) {
        ERROR("Release site " + rel_site_name + ": points_list item does not have three values.");
      }
      out << NAME_LOCATION << " = [" << point[0].asDouble() << ", " << point[1].asDouble() << ", " << point[2].asDouble() << "]";

      string orient = convert_orientation(release_site_item[KEY_ORIENT].asString());
      if (orient != "") {
        out << ", " << NAME_ORIENTATION << " = " << MDOT << NAME_ENUM_ORIENTATION << "." << orient;
      }

      out << ")";
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

  gen_param_id(out, NAME_REGION, region_expr, true);
}


static void error_release_pattern(const string& name) {
  ERROR("Requested release pattern " + name + " not found in data model.");
}


void PythonGenerator::generate_release_pattern(std::ostream& out, const std::string& name, std::string& delay_string) {

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

    // we found the release pattern
    check_version(KEY_RELEASE_PATTERN_LIST, release_pattern_item, VER_DM_2018_01_11_1330);

    gen_ctor_call(out, name, NAME_CLASS_RELEASE_PATTERN);
    gen_param(out, NAME_NAME, name, true);
    gen_param_expr(out, NAME_RELEASE_INTERVAL, release_pattern_item[KEY_RELEASE_INTERVAL], true);
    gen_param_expr(out, NAME_TRAIN_DURATION, release_pattern_item[KEY_TRAIN_DURATION], true);
    gen_param_expr(out, NAME_TRAIN_INTERVAL, release_pattern_item[KEY_TRAIN_INTERVAL], true);
    gen_param_expr(out, NAME_NUMBER_OF_TRAINS, release_pattern_item[KEY_NUMBER_OF_TRAINS], false);
    out << CTOR_END;

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

    gen_ctor_call(out, name, NAME_CLASS_RELEASE_SITE);
    gen_param(out, NAME_NAME, name, true);

    bool is_vol;
    if (shape != VALUE_LIST) {
      string species_name = release_site_item[KEY_MOLECULE].asString();
      if (!bng_mode) {
        gen_param_id(out, NAME_SPECIES, species_name, true);
      }
      else {
        gen_param(out, NAME_BNGL_SPECIES, remove_compartments(species_name), true);
      }


      string orientation = convert_orientation(release_site_item[KEY_ORIENT].asString());
      if (orientation != "") {
        // check that this is not a volume molecule
        is_vol = is_volume_species(species_name);
        if (is_vol) {
          cout <<
              "Ignoring orientation set for release site " << name << " with species " << species_name <<
              ", this species represent volume molecules.\n";
        }
        else {
          gen_param_enum(out, NAME_ORIENTATION, NAME_ENUM_ORIENTATION, orientation, true);
        }
      }
    }

    if (delay_string != "" && delay_string != "0") {
      gen_param_expr(out, NAME_RELEASE_TIME, delay_string, true);
    }

    if (rel_pat_name != "") {
      gen_param_id(out, NAME_RELEASE_PATTERN, rel_pat_name, true);
    }

    if (shape == VALUE_SPHERICAL) {
      gen_param_enum(out, NAME_SHAPE, NAME_ENUM_SHAPE, NAME_EV_SPHERICAL, true);
      gen_param_vec3(
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
      assert(molecule_list_name != "");
      bool diam_is_zero = release_site_item[KEY_SITE_DIAMETER] == "0";
      gen_param_id(out, NAME_MOLECULE_LIST, molecule_list_name, !diam_is_zero);
      if (!diam_is_zero) {
        gen_param_expr(out, NAME_SITE_DIAMETER, release_site_item[KEY_SITE_DIAMETER], false);
      }
    }
    else {
      ERROR("Shape " + shape + " is not supported yet");
    }

    if (shape != VALUE_LIST) {
      string quantity_type = release_site_item[KEY_QUANTITY_TYPE].asString();
      if (quantity_type == VALUE_NUMBER_TO_RELEASE) {
        gen_param_expr(out, NAME_NUMBER_TO_RELEASE, release_site_item[KEY_QUANTITY], false);
      }
      else if (quantity_type == VALUE_DENSITY) {
        string species_name = release_site_item[KEY_MOLECULE].asString();
        if (is_volume_species(species_name)) {
          gen_param_expr(out, NAME_CONCENTRATION, release_site_item[KEY_QUANTITY], false);
        }
        else {
          gen_param_expr(out, NAME_DENSITY, release_site_item[KEY_QUANTITY], false);
        }
      }
      else {
        ERROR("Quantity type " + quantity_type + " is not supported yet");
      }
    }

    out << CTOR_END;
    i++;
  }
}

} /* namespace MCell */
