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

/**
 * This file is directly included because it contains templates and also to avoid the
 * need to have two declarations (in .h + in .cpp) for each function.
 */

#ifndef SRC4_PYMCELLCONVERTER_GENERATOR_UTILS_H_
#define SRC4_PYMCELLCONVERTER_GENERATOR_UTILS_H_

#include <iostream>
#include <string>
#include <cassert>
#include <regex>

#include "libmcell/generated/gen_names.h"
#include "include/datamodel_defines.h"
#include "libmcell/api/api_utils.h"
#include "libmcell/api/python_export_constants.h"
#include "libmcell/api/python_export_utils.h"
#include "libmcell/api/api_common.h"

using namespace std;

namespace MCell {

using namespace API;
using Json::Value;

struct SharedGenData;

const uint TESTING_RXN_CLASS_CLEANUP_PERIODICITY = 100;
const uint TESTING_SPECIES_CLEANUP_PERIODICITY = 500;

const char* const NAME_PARAMETER = "parameter";

// using exception catching to recover from errors
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
      throw ConversionError(S("Error: Expected '") + #cond + "' is false. (" + __FUNCTION__ + " - " + __FILE__ + ":" + to_string(__LINE__) + ")"); \
    } \
  } while (0)

#define ERROR(msg) throw ConversionError(S("Error: ") + msg + " (function " + __FUNCTION__ + ")")


// key - MDL name, used in data model, value.first - Python name, value.second - BNGL name
const std::map<const std::string, std::pair<const std::string, const std::string>> mdl_functions_to_py_bngl_map {
  { "SQRT", {"math.sqrt", "sqrt"} },
  { "EXP", {"math.exp", "exp"} },
  { "LOG", {"math.log", "ln"} },
  { "LOG10", {"math.log10", "log10"} },
  { "SIN", {"math.sin", "sin"} },
  { "COS", {"math.cos", "cos"} },
  { "TAN", {"math.tan", "tan"} },
  { "ASIN", {"math.asin", "asin"} },
  { "ACOS", {"math.acos", "acos"} },
  { "ATAN", {"math.atan", "atan"} },
  { "ABS", {"abs", "abs"} },
  { "CEIL", {"math.ceil", "ceil"} },
  { "FLOOR", {"math.floor", "floor"} },
  { "MAX", {"max", "max"} },
  { "MIN", {"min", "min"} }
};


// when use_python_functions is true, function calls are replaced with Python function names
// when False, they are replaced with BNGL function names
std::string replace_function_calls_in_expr(const std::string& data_model_expr, const bool use_python_functions);

std::string get_module_name_w_prefix(const std::string& output_files_prefix, const std::string file_suffix);

void parse_rxn_rule_side(
    Json::Value& substances_node,
    std::vector<std::string>& substances,
    std::vector<std::string>& orientations);

// throws exception when the member is member is there
static Value& get_node(const string parent_name, Value& parent, const string name) {
  if (!parent.isMember(name)) {
    throw ConversionError("Error: Node '" + parent_name + "' does not contain expected node '" + name + "'.");
  }
  return parent[name];
}


// used when we know that the member is there
static Value& get_node(Value& parent, const string name) {
  assert(parent.isMember(name));
  return parent[name];
}


static string fix_dots_in_simple_species(const string& s) {
  string res = s;
  if (API::is_simple_species(s)) {
    replace(res.begin(), res.end(), '.', '_');
  }
  return res;
}


string remove_compartments(const std::string& species_name);

string get_single_compartment(const std::string& name, bool* has_multiple_compartments = nullptr);

string make_species_or_cplx(
    const SharedGenData& data,
    const std::string& name,
    const std::string& orient = "",
    const std::string& compartment = "");


static string make_species(const string bngl_str) {
  return S(MDOT) + API::NAME_CLASS_SPECIES + "('" + fix_dots_in_simple_species(bngl_str) + "')";
}


static void check_versions(
    const string node_name, Json::Value& node,
    const char* const version1, const char* const version2) {
  if (node[KEY_DATA_MODEL_VERSION].asString() != version1 &&
      node[KEY_DATA_MODEL_VERSION].asString() != version2) {
    throw ConversionError(
        "Error: version for " + node_name + " is " + node[KEY_DATA_MODEL_VERSION].asString() +
        ", expected " + version1 + " or " + version2 + ".");
  }
}


static void check_version(const string node_name, Json::Value& node, const char* const version) {
  if (!node.isMember(KEY_DATA_MODEL_VERSION) || node[KEY_DATA_MODEL_VERSION].asString() == "") {
    cout << "Warning: data model node " + node_name + " does not have its version specified, not checking the version.\n";
  }
  else if (node[KEY_DATA_MODEL_VERSION].asString() != version) {
    throw ConversionError(
        "Error: version for " + node_name + " is " + node[KEY_DATA_MODEL_VERSION].asString() +
        ", expected " + version + ".");
  }
}


static string convert_orientation(const string s, const bool return_any_orientation = false) {
  if (s == "\'") {
    return API::NAME_EV_UP;
  }
  else if (s == ",") {
    return API::NAME_EV_DOWN;
  }
  else if (s == ";" || s == "") {
    if (return_any_orientation && s == ";") {
      return API::NAME_EV_ANY;
    }
    else {
      return "";
    }
  }
  else {
    ERROR("Invalid orientation '" + s + "'.");
    return "INVALID_ORIENTATION";
  }
}

// NOTE: the same code is in mcell3_world_converter.cpp
static bool ends_with(std::string const & value, std::string const & ending)
{
    if (ending.size() > value.size()) {
      return false;
    }
    return std::equal(ending.rbegin(), ending.rend(), value.rbegin());
}

static string trim(const string& str)
{
    size_t first = str.find_first_not_of(' ');
    if (string::npos == first)
    {
        return str;
    }
    size_t last = str.find_last_not_of(' ');
    return str.substr(first, (last - first + 1));
}

string get_rxn_id(Json::Value& reaction_list_item, uint& unnamed_rxn_counter);

string create_count_name(
    const string& what_to_count, const string& where_to_count, const bool molecules_not_species);

uint get_num_counts_in_mdl_string(const string& mdl_string);

string remove_c_comment(const string& str);
string remove_whitespace(const string& str);

void process_single_count_term(
    const SharedGenData& data,
    const string& mdl_string,
    bool& rxn_not_mol,
    bool& molecules_not_species,
    string& what_to_count,
    string& where_to_count,
    string& orientation);

// sets val if the name_or_value is a floating point value,
// if not, tries to find the parameter and reads its value
// returns true on success
// parameters are not evaluated and only one level is tried,
// returns false if value was not obtained
bool get_parameter_value(Json::Value& mcell, const string& name_or_value, double& val);

bool is_volume_species(Json::Value& mcell, const std::string& species_name);

} // namespace MCell

#endif // SRC4_PYMCELLCONVERTER_GENERATOR_UTILS_H_
