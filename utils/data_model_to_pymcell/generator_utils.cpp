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
#include "api/python_export_constants.h"
#include "api/compartment_utils.h"

#include "generator_structs.h"

using namespace MCell::API;

namespace MCell {

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

  // cannot use BNGL parser directly because reactions may contain orientations

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
        if (isalnum(c) || c == '_' || c == '.' || c == '@') {
          state = CPLX;
          current_id = c;
        }
        else if (isblank(c)) {
          // ok
        }
        else {
          ERROR("Could not parse reaction side " + str + " (START).");
        }
        break;

      case CPLX:
        if (isalnum(c) || c == '_' || c == '.' || c == '@' || c == ':') {
          current_id += c;
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
        if (isalnum(c) || c == '_' || c == '@') {
          state = CPLX;
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
      if (c != ':') {
        res += c;
      }
    }
    else if (!in_compartment) {
      res += c;
    }
    i++;
  }
  return res;
}


string get_single_compartment(const std::string& name) {
  std::vector<std::string> compartments;
  API::get_compartment_names(name, compartments);

  if (compartments.empty()) {
    return "";
  }
  else {
    string res = compartments[0];
    for (size_t i = 1; i < compartments.size(); i++) {
      if (res != compartments[i]) {
        ERROR("Complexes may use a single compartment for now, error for " + name + ".");
      }
    }
    return res;
  }
}


static string make_cplx(const string bngl_str, const string orient = "", const string compartment = "") {
  string res = S(MDOT) + API::NAME_CLASS_COMPLEX + "('" + fix_dots_in_simple_species(bngl_str) + "'";
  if (orient != "") {
    res += S(", ") + API::NAME_ORIENTATION + " = " + MDOT + API::NAME_ENUM_ORIENTATION + "." + orient;
  }
  if (compartment != "") {
    res += S(", ") + API::NAME_COMPARTMENT_NAME + " = '" + compartment + "'";
  }
  res += ")";
  return res;
}


bool static is_mdot_superclass(const std::string& name) {
  return
      name == S(MDOT) + API::NAME_CV_AllMolecules ||
      name == S(MDOT) + API::NAME_CV_AllVolumeMolecules ||
      name == S(MDOT) + API::NAME_CV_AllSurfaceMolecules;
}

string make_species_or_cplx(
    const SharedGenData& data,
    const std::string& name,
    const std::string& orient,
    const std::string& compartment) {

  if (is_mdot_superclass(name)) {
    release_assert(orient == "");
    release_assert(compartment == "");
    return name;
  }

  if (!data.bng_mode) {
    stringstream ss;
    const SpeciesOrMolType* species_info = data.find_species_or_mol_type_info(name);

    if (species_info != nullptr && species_info->is_species) {
      // substance was declared as species, we can use its id directly
      ss << make_id(name) << "." << API::NAME_INST << "(";

      if (orient != "") {
        assert(compartment == "");
        ss <<
            API::NAME_ORIENTATION << " = " <<
            MDOT << API::NAME_ENUM_ORIENTATION << "." << orient;
      }

      if (compartment != "") {
        assert(orient == "");
        ss <<
            API::NAME_COMPARTMENT_NAME << " = '" << compartment << "'";
      }

      ss << ")";
      return ss.str();
    }
  }

  // otherwise we will generate a BNGL string
  return make_cplx(remove_compartments(name), orient, compartment);
}


string reaction_name_to_id(const string& json_name) {
  string res_name = json_name;
  replace(res_name.begin(), res_name.end(), ' ', '_');
  replace(res_name.begin(), res_name.end(), '.', '_');
  replace(res_name.begin(), res_name.end(), ')', '_');
  replace(res_name.begin(), res_name.end(), '(', '_');
  replace(res_name.begin(), res_name.end(), '!', '_');
  replace(res_name.begin(), res_name.end(), '~', '_');

  res_name = regex_replace(res_name, regex("<->"), "revto");
  res_name = regex_replace(res_name, regex("->"), "to");
  res_name = regex_replace(res_name, regex("\\+"), "plus");
  res_name = regex_replace(res_name, regex("\\?"), "any_bond");
  res_name = regex_replace(res_name, regex("'"), "_up");
  res_name = regex_replace(res_name, regex(","), "_down");
  res_name = regex_replace(res_name, regex(";"), "_any");
  res_name = regex_replace(res_name, regex("@"), "_at_");

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


string create_count_name(
    const string& what_to_count, const string& compartmnent, const string& where_to_count,
    const bool molecules_not_species) {

  // first remove all cterm_refixes
  regex pattern_cterm(COUNT_TERM_PREFIX);
  string what_to_count_no_cterm = regex_replace(what_to_count, pattern_cterm, "");
  string res = COUNT_PREFIX + fix_id(what_to_count_no_cterm);

  if (!molecules_not_species) {
    res += "_species";
  }

  if (compartmnent != "") {
    res += "_at_" + compartmnent;
  }

  if (where_to_count != WORLD && where_to_count != "") {
    res += "_" + where_to_count;
  }

  return res;
}


uint get_num_counts_in_mdl_string(const string& mdl_string) {
  uint res = 0;
  size_t pos = 0;
  while ((pos = mdl_string.find(COUNT, pos)) != string::npos) {
    res++;
    pos += strlen(COUNT);
  }
  return res;
}


string remove_c_comment(const string& str) {
  std::regex e ("/\\*.*\\*/");
  return std::regex_replace(str, e, "");
}

void process_single_count_term(
    const SharedGenData& data,
    const string& mdl_string,
    bool& rxn_not_mol,
    bool& molecules_not_species,
    string& what_to_count,
    string& compartment,
    string& where_to_count,
    string& orientation) {

  // mdl_string is always in the form COUNT[what,where]
  size_t start_brace = mdl_string.find('[');
  size_t comma = mdl_string.rfind(',');
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

  size_t pos_at = what_to_count.find('@');
  compartment = ""; // FIXME: use get_single_compartment & remove_compartment
  if (pos_at != string::npos) {
    compartment = what_to_count.substr(pos_at + 1);
    what_to_count = what_to_count.substr(0, pos_at);
  }

  // default is 'molecules_pattern', for now we are storing the
  // as a comment because the counting type belongs to the count term,
  // not to the whole 'reaction_output_list'
  // TODO: this should resolved in a better way
  molecules_not_species = true;
  if (what_to_count.find(MARKER_SPECIES_COMMENT) != string::npos) {
    molecules_not_species = false;
  }

  what_to_count = remove_c_comment(what_to_count);

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

  if (data.find_reaction_rule_info(what_to_count) != nullptr) {
    rxn_not_mol = true;
  }
  else if (data.bng_mode ||
      data.find_species_or_mol_type_info(what_to_count) != nullptr ||
      !API::is_simple_species(what_to_count)) {
    // if we did not find the name to be a reaction, we assume it is simple species or complex pattern
    rxn_not_mol = false;
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
  // where_to_count can now look like this: "Cube[ALL"
  size_t brace = where_to_count.find('[');
  if (brace != string::npos) {
    if (where_to_count.substr(brace + 1) == REGION_ALL_NAME) {
      where_to_count = where_to_count.substr(0, brace);
    }
    else {
      where_to_count[brace] = '_';
    }
  }

  if (compartment != "" && where_to_count != "") {
      ERROR("Cannot both specify location and compartment for a count location '" + mdl_string + "'.");
  }
}

// sets val if the name_or_value is a floating point value,
// if not, tries to find the parameter and reads its value
// returns true on success
// parameters are not evaluated and only one level is tried,
// returns false if value was not obtained
bool get_parameter_value(Json::Value& mcell, const string& name_or_value, double& val) {
  try {
    val = stod(name_or_value);
    return true;
  }
  catch (const std::invalid_argument& ia) {
    // not a float, try to get parameter value
    if (mcell.isMember(KEY_PARAMETER_SYSTEM) && mcell.isMember(KEY_MODEL_PARAMETERS)) {
      Json::Value& params = mcell[KEY_PARAMETER_SYSTEM][KEY_MODEL_PARAMETERS];
      for (Value::ArrayIndex i = 0; i < params.size(); i++) {
        if (params[i][KEY_NAME].asString() == name_or_value) {
          try {
            if (params[i].isMember(KEY__EXTRAS)) {
              val = stod(params[i][KEY__EXTRAS][KEY_PAR_VALUE].asString());
            }
            else {
              val = stod(params[i][KEY_PAR_EXPRESSION].asString());
            }
            return true;
          }
          catch (const std::invalid_argument& ia) {
            return false;
          }
        }
      }
    }
  }
  return false;
}


} // namespace MCell
