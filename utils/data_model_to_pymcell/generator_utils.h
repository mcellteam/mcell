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

/**
 * This file is directly included because it contains templates and also to avoid the
 * need to have two declarations (in .h + in .cpp) for each function.
 */

#ifndef SRC4_PYMCELLCONVERTER_GENERATOR_UTILS_H_
#define SRC4_PYMCELLCONVERTER_GENERATOR_UTILS_H_

#include <exception>
#include <iostream>
#include <string>
#include <cassert>

#include "../../libmcell/generated/gen_names.h"
#include "../../include/datamodel_defines.h"

#include "generator_constants.h"

using namespace std;

namespace MCell {

using Json::Value;

typedef invalid_argument ConversionError;

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

#define ERROR(msg) throw ConversionError(msg)


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


static void gen_comma(ostream& out, Value::ArrayIndex index, Value& array) {
  if (index + 1 != array.size()) {
    out << ", ";
  }
}


template<typename T>
static void gen_comma(ostream& out, size_t index, const vector<T>& array) {
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


static string make_enum_value(const string enum_name, const string value) {
  return MDOT + enum_name + "." + value;
}


static void check_version(const string node_name, Json::Value& node, const char* const version) {
  if (node[KEY_DATA_MODEL_VERSION].asString() != version) {
    throw ConversionError(
        "Error: version for " + node_name + " is " + node[KEY_DATA_MODEL_VERSION].asString() +
        ", expected" + version);
  }
}


void gen_ctor_call(ofstream& out, string name, string class_name, bool has_params = true) {
  out << name << " = " << MDOT << class_name;
  if (has_params) {
    out << "(\n";
  }
  else {
    out << "()\n";
  }
}


void gen_method_call(ofstream& out, string obj, string method, string param = "") {
  out << obj << "." << method << "(" << param << ")\n";
}


template<typename T>
void gen_param(ofstream& out, string name, T value, bool comma) {
  out << IND << name << " = " << value << (comma?",":"") << "\n";
}


template<>
void gen_param(ofstream& out, string name, Json::Value& value, bool comma) {
  out << IND << name << " = '" << value.asString() << "'" << (comma?",":"") << "\n";
}


template<>
void gen_param(ofstream& out, string name, string value, bool comma) {
  out << IND << name << " = '" << value << "'" << (comma?",":"") << "\n";
}


template<>
void gen_param(ofstream& out, string name, int value, bool comma) {
  out << IND << name << " = " << value << (comma?",":"") << "\n";
}


template<>
void gen_param(ofstream& out, string name, double value, bool comma) {
  out << IND << name << " = " << value << (comma?",":"") << "\n";
}


void gen_param_id(ofstream& out, string name, string id, bool comma) {
  out << IND << name << " = " << id << (comma?",":"") << "\n";
}


void gen_param_id(ofstream& out, string name, Json::Value& id, bool comma) {
  out << IND << name << " = " << id.asString() << (comma?",":"") << "\n";
}

void gen_param_int(ofstream& out, string name, Json::Value& value, bool comma) {
  string s = value.asString();
  int v = stoi(s);
  gen_param(out, name, v, comma);
}


void gen_param_double(ofstream& out, string name, Json::Value& value, bool comma) {
  string s = value.asString();
  double v = stod(s);
  gen_param(out, name, v, comma);
}


void gen_param_enum(ofstream& out, string name, string enum_name, string enum_value, bool comma) {
  out << IND << name << " = " << make_enum_value(enum_name, enum_value) << (comma?",":"") << "\n";
}

void gen_param_list(ofstream& out, string name, const vector<string>& values, bool comma) {
  out << IND << name << " = [";
  for (size_t i = 0; i < values.size(); i++) {
    out << values[i];
    gen_comma(out, i , values);
  }

  out << "]" << (comma?",":"") << "\n";
}

void gen_param_vec3(ofstream& out, string name, Json::Value& x, Json::Value& y, Json::Value& z, bool comma) {
  out << IND <<
      name << " = " << MDOT << VEC3 <<
      "(" << x.asString() << ", " << y.asString() << ", " << z.asString() << ")" << (comma?",":"") << "\n";
}


template<typename T>
void gen_assign(ofstream& out, string obj_name, string field_name1, string field_name2, T value) {
  out << obj_name << "." << field_name1 << "." << field_name2 << " = " << value << "\n";
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
    out << substances[i] << "." << API::NAME_INST << "()";
    print_comma(out, i, substances);
  }
  out << " ]";
}

} // namespace MCell

#endif // SRC4_PYMCELLCONVERTER_GENERATOR_UTILS_H_
