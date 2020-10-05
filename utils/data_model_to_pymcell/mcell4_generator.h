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

// TODO: add assert after every find_* call or add a call that checks that
//       the thing we are searching for exists

#ifndef SRC4_MCELL4_GENERATOR_H_
#define SRC4_MCELL4_GENERATOR_H_

#include <string>
#include <fstream>
#include "json/json.h"

#include "generator_structs.h"
#include "python_generator.h"
#include "bngl_generator.h"

// TODO: use a different namespace, we are not using anything from MCell
// except for constants
namespace MCell {

class MCell4Generator {
public:
  bool generate(
      const std::string& input_file,
      const std::string& output_files_prefix_,
      const bool bng_mode_,
      const bool debug_mode_,
      const bool cellblender_viz
  );

private:
  void reset();

  std::string get_module_name(const std::string file_suffix);
  std::string make_import(const std::string file_suffix);

  void open_and_check_file(
      const std::string file_suffix, std::ofstream& out,
      const bool for_append = false,
      const bool bngl = false);

  void check_scripting();

  void generate_parameters();

  void generate_species_and_mol_types(std::ostream& out, std::vector<SpeciesOrMolType>& species_and_mt_info);

  void generate_variable_rate(const std::string& rate_array_name, Json::Value& variable_rate_text);
  std::vector<IdLoc> generate_reaction_rules(std::ostream& out);

  void generate_subsystem();

  std::string generate_single_geometry_object(
      std::ofstream& out, const int index, Json::Value& object);
  std::vector<std::string> generate_geometry();

  void generate_surface_classes_assignment(std::ostream& out);
  void generate_instantiation(const std::vector<std::string>& geometry_objects);

  void generate_observables(const bool cellblender_viz);

  void generate_config(std::ostream& out);
  void generate_model(const bool print_failed_marker);

private:
  BNGLGenerator* bng_gen;
  std::ofstream bng_out;

  PythonGenerator* python_gen;

  // parameters, subsystem, and instantiation are always generated
  bool geometry_generated;
  bool observables_generated;

  SharedGenData data;
};

} /* namespace MCell */

#endif /* SRC4_MCELL4_GENERATOR_H_ */
