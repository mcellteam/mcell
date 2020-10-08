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
#include "generator_constants.h"

namespace MCell {

string get_filename(const string& output_files_prefix, const string file_suffix, const char* ext) {
  if (output_files_prefix == "" || output_files_prefix.back() == '/' || output_files_prefix.back() == '\\') {
    return output_files_prefix + file_suffix + ext;
  }
  else {
    return output_files_prefix + "_" + file_suffix + ext;
  }
}


void open_and_check_file_w_prefix(
    const string& output_files_prefix, const string file_suffix, ofstream& out,
    const bool for_append, const bool bngl) {

  string file_name = get_filename(output_files_prefix, file_suffix, (bngl) ? BNGL_EXT : PY_EXT);

  if (for_append) {
    cout << "Appending to " + file_name + ".\n";
    out.open(file_name, ios::app);
  }
  else {
    cout << "Generating file " + file_name + ".\n";
    out.open(file_name);
  }
  out.precision(FLOAT_OUT_PRECISION);
  if (!out.is_open()) {
    throw ConversionError("Could not open file '" + file_name + "' for writing.");
  }
}


std::string get_module_name_w_prefix(const std::string& output_files_prefix, const std::string file_suffix) {

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


void parse_rxn_rule_side(
    Json::Value& substances_node,
    std::vector<std::string>& substances,
    std::vector<std::string>& orientations) {

  // cannot use BNGL parser because reactions may contain orientations

  substances.clear();
  orientations.clear();

  string str = substances_node.asString();

  // finite automata to parse the reaction side string, e.g. "a + b"
  enum state_t {
    START,
    CPLX,
    IN_PAREN,
    AFTER_MOL,
    AFTER_ORIENT,
    AFTER_PLUS
  };

  state_t state = START;
  string current_id;
  for (size_t i = 0; i < str.size(); i++) {
    char c = str[i];
    switch (state) {
      case START:
        if (isalnum(c) || c == '_') {
          state = CPLX;
          current_id = c;
        }
        else if (c == '.') {
          state = CPLX;
          current_id = '_';
        }
        else if (isblank(c)) {
          // ok
        }
        else {
          ERROR("Could not parse reaction side " + str + " (START).");
        }
        break;

      case CPLX:
        if (isalnum(c) || c == '_' || c == '.') {
          current_id += c;
        }
        else if (c == '.') {
          state = CPLX;
          current_id += '_';
        }
        else if (c == '(') {
          state = IN_PAREN;
          current_id += '(';
        }
        else if (isblank(c) || c == '+' || c == '\'' || c == ',' || c == ';') {
          substances.push_back(current_id);
          orientations.push_back("");
          if (c == '\'' || c == ',' || c == ';') {
            orientations.back() = c;
          }
          current_id = "";
          if (c == '+') {
            state = AFTER_PLUS;
          }
          else {
            state = AFTER_MOL;
          }
        }
        else {
          ERROR("Could not parse reaction side " + str + " (MOL_ID).");
        }
        break;

      case IN_PAREN:
        if (isalnum(c) || c == '_' || c == ',' || c == '~' || c == '!' || c == '+' || c == '?') {
          current_id += c;
        }
        else if (c == ')') {
          state = CPLX;
          current_id += ')';
        }
        else {
          ERROR("Could not parse reaction side " + str + " (IN_PAREN).");
        }
        break;

      case AFTER_MOL:
        if (c == '+') {
          state = AFTER_PLUS;
        }
        else if (c == '\'') {
          state = AFTER_ORIENT;
          orientations.back() = c;
        }
        else if (c == ',') {
          state = AFTER_ORIENT;
          orientations.back() = c;
        }
        else if (c == ';') {
          state = AFTER_ORIENT;
          orientations.back() = c;
        }
        else if (isblank(c)) {
          // ok
        }
        else {
          ERROR("Could not parse reaction side " + str + " (AFTER_MOL).");
        }
        break;

      case AFTER_ORIENT:
        if (c == '+') {
          state = AFTER_PLUS;
        }
        else if (isblank(c)) {
          // ok
        }
        else {
          ERROR("Could not parse reaction side " + str + " (AFTER_ORIENT).");
        }
        break;

      case AFTER_PLUS:
        if (isalnum(c) || c == '_') {
          state = CPLX;
          current_id = c;
        }
        else if (c == '.') {
          state = CPLX;
          current_id = '_';
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
    orientations.push_back("");
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


static string make_cplx_inst(const string bngl_str, const string orient = "") {
  string res = S(MDOT) + API::NAME_CLASS_COMPLEX_INSTANCE + "('" + fix_dots_in_simple_species(bngl_str) + "'";
  if (orient != "") {
    res += S(", ") + API::NAME_ORIENTATION + " = " + MDOT + API::NAME_ENUM_ORIENTATION + "." + orient;
  }
  res += ")";
  return res;
}


string make_species_or_cplx_inst(
    const SharedGenData& data, const std::string& name, const std::string& orient) {

  if (!data.bng_mode) {
    stringstream ss;
    const SpeciesOrMolType* species_info = data.find_species_or_mol_type_info(name);

    if (species_info != nullptr && species_info->is_species) {
      // substance was declared as species, we can use its id directly
      ss << make_id(name) << "." << API::NAME_INST << "(";

      if (orient != "") {
        ss <<
            API::NAME_ORIENTATION << " = " <<
            MDOT << API::NAME_ENUM_ORIENTATION << "." << orient;
      }

      ss << ")";
      return ss.str();
    }
  }

  // otherwise we will generate a BNGL string
  return make_cplx_inst(remove_compartments(name), orient);
}


string fix_id(const std::string& str) {
  string res;
  for (char c: str) {
    if (c == '+') {
      res += "_plus_";
    }
    else if (c == '-') {
      res += "_minus_";
    }
    else if (c == '?') {
      res += "_anybond_";
    }
    else if (c == '!') {
      res += "_bond_";
    }
    else if (c == '(') {
      res += "_ps_";
    }
    else if (c == ')') {
      res += "_pe_";
    }
    else if (
        c == ' ' || c == '.' || c == '_' ||
        c == ',' || c == '~') {
      res += "_";
    }
    else if (isalnum(c)) {
      res += c;
    }
    // ignoring the rest of the characters
  }
  return res;
}



string reaction_name_to_id(const string& json_name) {
  string res_name = json_name;
  replace(res_name.begin(), res_name.end(), ' ', '_');
  replace(res_name.begin(), res_name.end(), '.', '_');
  replace(res_name.begin(), res_name.end(), ')', '_');
  replace(res_name.begin(), res_name.end(), '(', '_');
  replace(res_name.begin(), res_name.end(), '!', '_');

  res_name = regex_replace(res_name, regex("<->"), "revto");
  res_name = regex_replace(res_name, regex("->"), "to");
  res_name = regex_replace(res_name, regex("\\+"), "plus");
  res_name = regex_replace(res_name, regex("'"), "_up");
  res_name = regex_replace(res_name, regex(","), "_down");
  res_name = regex_replace(res_name, regex(";"), "_any");

  return res_name;
}


string get_rxn_id(Json::Value& reaction_list_item, uint& unnamed_rxn_counter) {
  string name = reaction_list_item[KEY_RXN_NAME].asString();
  if (name == "") {
    name = UNNAMED_REACTION_RULE_PREFIX + to_string(unnamed_rxn_counter);
    unnamed_rxn_counter++;
  }
  else {
    name = reaction_name_to_id(name);
  }
  return name;
}


} // namespace MCell
