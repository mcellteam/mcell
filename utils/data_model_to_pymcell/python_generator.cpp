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

  sc_names.clear();

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

} /* namespace MCell */
