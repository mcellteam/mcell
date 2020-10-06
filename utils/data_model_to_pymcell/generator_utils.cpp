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

  substances.clear();
  orientations.clear();

  string str = substances_node.asString();

  // finite automata to parse the reaction side string, e.g. "a + b"
  enum state_t {
    START,
    ID,
    AFTER_ID,
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
          state = ID;
          current_id = c;
        }
        else if (c == '.') {
          state = ID;
          current_id = '_';
        }
        else if (isblank(c)) {
          // ok
        }
        else {
          ERROR("Could not parse reaction side " + str + " (START).");
        }
        break;

      case ID:
        if (isalnum(c) || c == '_') {
          current_id += c;
        }
        else if (c == '.') {
          state = ID;
          current_id += '_';
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
            state = AFTER_ID;
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
          ERROR("Could not parse reaction side " + str + " (AFTER_ID).");
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
          ERROR("Could not parse reaction side " + str + " (AFTER_ID).");
        }
        break;

      case AFTER_PLUS:
        if (isalnum(c) || c == '_') {
          state = ID;
          current_id = c;
        }
        else if (c == '.') {
          state = ID;
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

} // namespace MCell
