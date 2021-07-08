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

#ifndef UTILS_DATA_MODEL_TO_PYMCELL_PYTHON_GENERATOR_H_
#define UTILS_DATA_MODEL_TO_PYMCELL_PYTHON_GENERATOR_H_

#include <iostream>
#include <string>
#include <vector>

#include "generator_structs.h"

namespace MCell {

class MCell4Generator;

class PythonGenerator {
public:
  PythonGenerator(SharedGenData& data_)
    : data(data_), mcell(data_.mcell), unnamed_surf_class_counter(0)
    {
  }

  void generate_parameters(std::ostream& out);

  void generate_species_and_mol_types(std::ostream& out, std::vector<SpeciesOrMolType>& species_and_mt_info);
  void generate_surface_classes(std::ostream& out, std::vector<std::string>& sc_names);

  // the parameters file must be closed because we might append some code to it
  std::string generate_single_reaction_rule(std::ostream& out, Json::Value& reaction_list_item);
  void generate_reaction_rules(std::ostream& out, std::vector<IdLoc>& rxn_names);

  void generate_geometry(std::ostream& out, std::vector<std::string>& geometry_objects);

  void generate_release_sites(std::ostream& out, std::vector<std::string>& release_site_names);
  void generate_surface_classes_assignments(std::ostream& out);
  void generate_compartment_assignments(std::ostream& out);

  void generate_viz_outputs(std::ostream& out, const bool cellblender_viz, std::vector<std::string>& viz_output_names);

  void generate_all_bngl_reaction_rules_used_in_observables(std::ostream& out);

  void generate_single_count(
      std::ostream& out,
      const std::string& count_name,
      const std::string& observable_name,
      const std::string& file_name,
      const std::string& count_term_name,
      const std::string& mul_div_str,
      const std::string& rxn_step);

  std::string generate_count_terms_for_expression(
        ostream& out,
        const string& mdl_string, // may be empty, in that case we use what_to_count and where_to_count
        const std::string& what_to_count,
        const std::string& where_to_count,
        const std::string& orientation,
        const bool rxn_not_mol);

  std::string generate_single_molecule_release_info_array(
      std::ostream& out,
      std::string& rel_site_name,
      Json::Value& release_site_list,
      Json::Value::ArrayIndex begin,
      Json::Value::ArrayIndex end);

  void generate_release_pattern(std::ostream& out, const std::string& name, std::string& delay_string);

private:
  void generate_single_parameter(std::ostream& out, Json::Value& parameter);

  std::string generate_component_type(
      std::ostream& out, Json::Value& bngl_component_item, const std::string& mol_type_name);

  std::string generate_single_species_or_mol_type(
      std::ostream& out, Json::Value& molecule_list_item,
      const bool generate_species, const std::vector<std::string>& component_names = std::vector<std::string>());

  SpeciesOrMolType generate_single_species_or_mol_type_w_components(
      std::ostream& out, Json::Value& molecule_list_item);


  std::string generate_single_count_term(
      ostream& out,
      const std::string& what_to_count,
      const std::string& where_to_count,
      const std::string& orientation,
      const bool molecules_not_species,
      const bool rxn_not_mol);

  void generate_rxn_rule_side(std::ostream& out, Json::Value& substances_node);

  void get_surface_class_property_info(
      const string& sc_name, Json::Value& property,
      std::string& name, std::string& type_name,
      std::string& affected_mols, std::string& orientation, std::string& clamp_concentration);

  void generate_variable_rate(const std::string& rate_array_name, Json::Value& variable_rate_text);

  std::string generate_single_geometry_object(
      std::ostream& out, const int index, Json::Value& object);

  bool is_volume_mol_type(const std::string& mol_type_name);

  std::vector<std::string> get_species_to_visualize();


private:
  SharedGenData& data; // owned by MCell4Generator
  Json::Value& mcell;
  uint unnamed_surf_class_counter;
};

} /* namespace MCell */

#endif /* UTILS_DATA_MODEL_TO_PYMCELL_PYTHON_GENERATOR_H_ */
