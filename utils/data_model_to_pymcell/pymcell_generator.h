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

#ifndef SRC4_PYMCELLCONVERTER_H_
#define SRC4_PYMCELLCONVERTER_H_

#include <string>
#include <fstream>

#include "json/json.h"

// use different namespace?
namespace MCell {

class PymcellGenerator {
public:
  bool generate(
      const std::string& input_file,
      const std::string& output_file_prefix_
  );

private:
  void reset();

  std::string get_filename(const std::string file_suffix);
  std::string get_module_name(const std::string file_suffix);
  std::string make_import(const std::string file_suffix);

  void open_and_check_file(const std::string file_suffix, std::ofstream& out);

  // TODO: shorten the name to gen?
  void generate_parameters();

  std::vector<std::string> generate_species(std::ofstream& out);
  std::vector<std::string> generate_reaction_rules(std::ofstream& out);
  void generate_subsystem();

  std::string generate_single_geometry_object(
      std::ofstream& out, const int index, Json::Value& object);
  std::vector<std::string> generate_geometry();

  std::vector<std::string> generate_release_sites(std::ofstream& out);
  void generate_instantiation(const std::vector<std::string>& geometry_objects);

  std::vector<std::string> generate_viz_outputs(std::ofstream& out);
  std::vector<std::string> generate_counts(std::ofstream& out);
  void generate_observables();

  void generate_config(std::ofstream& out);
  void generate_model();


  std::string output_files_prefix;

  // parameters, subsystem, and instantiation are always generated
  bool geometry_generated;
  bool observables_generated;

  uint unnamed_rxn_counter;

  std::vector<std::string> all_species_names;

  // mcell node of the loaded JSON file
  Json::Value mcell;
};

} /* namespace MCell */

#endif /* SRC4_PYMCELLCONVERTER_H_ */