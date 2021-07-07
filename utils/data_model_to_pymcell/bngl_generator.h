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

#ifndef UTILS_DATA_MODEL_TO_PYMCELL_BNGL_GENERATOR_H_
#define UTILS_DATA_MODEL_TO_PYMCELL_BNGL_GENERATOR_H_

#include <string>
#include <iostream>

#include "json/json.h"
#include "bng/bngl_names.h"

namespace MCell {

typedef std::map<std::string, std::set<std::string>> ParentToChildCompartmentsMap;

class BNGLGenerator {
public:
  BNGLGenerator(
      const std::string& bngl_filename_,
      std::ostream& bng_out_,
      SharedGenData& data_)
    : bngl_filename(bngl_filename_), bng_out(bng_out_), data(data_) {
  }

  void generate_units_information_header();
  void generate_parameters(std::ostream& python_out);
  void generate_mol_types(std::ostream& python_out);

  void generate_compartments();

  void open_reaction_rules_section() { bng_out << BNG::BEGIN_REACTION_RULES << "\n"; }
  std::string generate_single_reaction_rule(Json::Value& reaction_list_item, const bool generate_name);
  void close_reaction_rules_section() { bng_out << BNG::END_REACTION_RULES << "\n\n"; }

  bool can_express_count_with_bngl(
      const bool single_term,
      const bool rxn_not_mol,
      const std::string& where_to_count,
      const std::string& orientation,
      const std::string& mul_div_str,
      const std::string& rxn_step) const;

  void open_observables_section() { bng_out << BNG::BEGIN_OBSERVABLES << "\n"; }
  void generate_single_count(
      const std::string& observable_name,
      const std::string& what_to_count,
      const bool molecules_not_species);
  void close_observables_section() { bng_out << BNG::END_OBSERVABLES << "\n\n"; }

  bool can_express_release_with_bngl(Json::Value& release_site_item);

  void open_seed_species_section() { bng_out << BNG::BEGIN_SEED_SPECIES << "\n"; }
  void generate_single_release_site(
      const std::string& bngl_cplx,
      const std::string& quantity);
  void close_seed_species_section() { bng_out << BNG::END_SEED_SPECIES << "\n\n"; }


  void add_comment(const std::string& text) { bng_out << "# " << text << "\n"; }

private:
  void generate_single_bngl_parameter(Json::Value& parameter);
  void generate_single_python_parameter(std::ostream& python_out, Json::Value& parameter);

  void generate_bngl_mol_type(Json::Value& molecule_list_item);
  void generate_python_mol_type_info(std::ostream& python_out, Json::Value& molecule_list_item);

  void get_compartment_volume_and_area(const std::string& name, double& volume, double& area);
  Json::Value& find_geometry_object(const std::string& name);
  void generate_single_compartment(Json::Value& model_object);

  const std::string bngl_filename;
  std::ostream& bng_out;
  SharedGenData& data; // owned by MCell4Generator

  std::map<std::string, std::pair<double, double>> compartment_name_to_volume_area_cache;

  ParentToChildCompartmentsMap volume_compartment_children;
};

} /* namespace MCell */

#endif /* UTILS_DATA_MODEL_TO_PYMCELL_BNGL_GENERATOR_H_ */
