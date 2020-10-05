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

#ifndef UTILS_DATA_MODEL_TO_PYMCELL_BNGL_GENERATOR_H_
#define UTILS_DATA_MODEL_TO_PYMCELL_BNGL_GENERATOR_H_

#include <string>
#include <iostream>

#include "json/json.h"

namespace MCell {

class BNGLGenerator {
public:
  BNGLGenerator(
      const std::string& bngl_filename_,
      std::ostream& bng_out_,
      Json::Value& mcell_,
      const std::string& output_files_prefix_,
      uint& unnamed_rxn_counter_)
    : bngl_filename(bngl_filename_), bng_out(bng_out_),
      mcell(mcell_), output_files_prefix(output_files_prefix_), unnamed_rxn_counter(unnamed_rxn_counter_) {
  }

  void generate_parameters(std::ostream& python_out);
  void generate_mol_types(std::ostream& python_out);

  void open_reaction_rules_section() { bng_out << "begin reaction rules\n"; }
  std::string generate_single_reaction_rule(Json::Value& reaction_list_item, const bool generate_name);
  void close_reaction_rules_section() { bng_out << "end reaction rules\n"; }
  void generate_python_decl_bngl_rxn_rule(std::ostream& python_out, const std::string& name);

  void add_comment(const std::string& text) { bng_out << "# " << text << "\n"; }

private:
  void generate_single_bngl_parameter(Json::Value& parameter);
  void generate_single_python_parameter(std::ostream& python_out, Json::Value& parameter);

  void generate_bngl_mol_type(Json::Value& molecule_list_item);
  void generate_python_mol_type_info(std::ostream& python_out, Json::Value& molecule_list_item);

  const std::string& bngl_filename;
  std::ostream& bng_out;
  Json::Value& mcell;
  const std::string& output_files_prefix;
  uint& unnamed_rxn_counter; // owned by MCell4Generator
};

} /* namespace MCell */

#endif /* UTILS_DATA_MODEL_TO_PYMCELL_BNGL_GENERATOR_H_ */
