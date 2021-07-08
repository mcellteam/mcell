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
  bool generate(const SharedGenData& opts);

private:
  void reset();

  std::string get_module_name(const std::string file_suffix);
  std::string make_import(const std::string file_suffix);

  void open_and_check_file(
      const std::string file_suffix, std::ofstream& out,
      const bool for_append = false,
      const bool bngl = false);

  void generate_customization();

  void generate_shared();

  void generate_simulation_setup_parameter(std::ostream& out, const string& name, const string& value);
  void generate_parameters();

  std::string generate_species_and_mol_types(std::ostream& out, std::vector<SpeciesOrMolType>& species_and_mt_info);

  void generate_variable_rate(const std::string& rate_array_name, Json::Value& variable_rate_text);
  std::vector<IdLoc> generate_reaction_rules(std::ostream& out, const std::vector<std::string>& surf_class_names);

  void find_required_compartments(std::set<std::string>& compartments);
  void analyze_and_generate_bngl_compartments(std::ostream& out);
  void generate_subsystem();

  std::string generate_single_geometry_object(
      std::ofstream& out, const int index, Json::Value& object);
  std::vector<std::string> generate_geometry();

  void generate_release_sites(std::ostream& out, std::vector<std::string>& release_site_names);
  void generate_instantiation(const std::vector<std::string>& geometry_objects);

  void generate_counts(std::ostream& out, std::vector<std::string>& python_counts, bool& has_bng_observables);
  void generate_observables();

  void generate_config(std::ostream& out);
  void generate_model(const bool print_failed_marker);
  void generate_customization_template();
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
