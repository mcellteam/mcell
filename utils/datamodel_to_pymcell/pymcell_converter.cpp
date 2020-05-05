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

#include <exception>
#include <iostream>
#include <string>
#include <cassert>

#include "pymcell_converter.h"

#include "../../libmcell/generated/gen_names.h"
#include "../../include/datamodel_defines.h"

using namespace std;
using namespace MCell::API;

namespace MCell {

using Json::Value;

typedef std::invalid_argument ConversionError;


const char* const GEOMETRY_SUFFIX = "_geometry";
const char* const PY_EXT = ".py";

const char* const MDOT = "m.";
const char* const MCELL_IMPORT = "import mcell as m\n\n";

const char* const IND4 = "    ";
const char* const BLOCK_BEGIN1 = "# ---- ";
const char* const BLOCK_BEGIN2 = " ----\n";
const char* const BLOCK_END1 = "# ^^^^ ";
const char* const BLOCK_END2 = " ^^^^\n\n";

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
      throw ConversionError("Expected '" + #cond + "' is false. (" + __FUNCTION__ + " - " + __FILE__ + ":" + __LINE__ + ")"); \
    } \
  } while (0)


// auxiliary method to simply convert to std::string for when concatenanting string
static string S(const char* s) {
  return string(s);
}

static void check_file(ofstream& out, const std::string& file_name) {
  if (!out.is_open()) {
    throw ConversionError("Could not open file '" + file_name + "' for writing.");
  }
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
static void print_comma(ostream& out, size_t index, const vector<T>& array) {
  if (index + 1 != array.size()) {
    out << ", ";
  }
}

// the aim is to generate as much of output as possible,
// therefore we are using exceptions
bool PymcellConverter::convert(
    const string& input_file,
    const string& output_files_prefix_
) {
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

  CHECK(convert_geometry(), failed);

  return !failed;
}


void PymcellConverter::convert_single_geomerty_object(
    ostream& out, const int index, Value& object) {

  string parent_name = S(KEY_OBJECT_LIST) + "[" + to_string(index) + "]";

  string name = get_node(parent_name, object, KEY_NAME).asString();
  Value& vertex_list = get_node(parent_name, object, KEY_VERTEX_LIST);
  Value& element_connections = get_node(parent_name, object, KEY_ELEMENT_CONNECTIONS);
  // TODO: material_names

  out << BLOCK_BEGIN1 << name << BLOCK_BEGIN2;

  // vertex_list
  string id_vertex_list = name + "_" + NAME_VERTEX_LIST;
  out << id_vertex_list << " = [\n";
  for (Value::ArrayIndex i = 0; i < vertex_list.size(); i++) {
    out << IND4 << "[";
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
    out << IND4 << "[";
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
          out << "\n" << IND4;
        }
        out << include_elements[k].asInt();
        print_comma(out, k, include_elements);
      }
      out << "\n] #" << sr_element_connections_name << "\n\n";

      out << sr_global_name << " = " << MDOT << NAME_CLASS_SURFACE_REGION << "(\n";
      out << IND4 << NAME_NAME << " = '" << name << "',\n";
      out << IND4 << NAME_ELEMENT_CONNECTIONS << " = " << sr_element_connections_name << "\n";
      out << ")\n\n";

      sr_global_names.push_back(sr_global_name);
    }
  }

  // object creation itself
  out << name << " = " << MDOT << NAME_CLASS_GEOMETRY_OBJECT << "(\n";
  out << IND4 << NAME_NAME << " = '" << name << "',\n";
  out << IND4 << NAME_VERTEX_LIST << " = " << id_vertex_list << ",\n";
  out << IND4 << NAME_ELEMENT_CONNECTIONS << " = " << id_element_connections << ",\n";
  out << IND4 << NAME_SURFACE_REGIONS << " = [";
  for (size_t i = 0; i < sr_global_names.size(); i++) {
    out << sr_global_names[i];
    print_comma(out, i, sr_global_names);
  }
  out << "]\n)\n";

  out << BLOCK_END1 << name << BLOCK_END2;
}



void PymcellConverter::convert_geometry() {

  if (!mcell.isMember(KEY_GEOMETRICAL_OBJECTS)) {
    return;
  }
  Value& geometrical_objects = get_node(mcell, KEY_GEOMETRICAL_OBJECTS);
  if (!geometrical_objects.isMember(KEY_OBJECT_LIST)) {
    return;
  }
  Value& object_list = get_node(geometrical_objects, KEY_OBJECT_LIST);

  string file_name = output_files_prefix + GEOMETRY_SUFFIX + PY_EXT;
  ofstream out;
  out.open(file_name);
  check_file(out, file_name);

  out << MCELL_IMPORT;

  for (Value::ArrayIndex i = 0; i < object_list.size(); i++) {
    Value& object = object_list[i];
    convert_single_geomerty_object(out, i, object);
  }

  out.close();
  cout << "Created file " + file_name + ".\n";
}


} /* namespace MCell */
