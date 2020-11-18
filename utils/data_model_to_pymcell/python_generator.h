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
      const std::string& name,
      const std::string& mdl_file_prefix_orig,
      const std::string& what_to_count,
      const std::string& compartment,
      const std::string& where_to_count, // empty for WORLD
      const std::string& orientation,
      const std::string& multiplier_str,
      const std::string& rxn_step,
      const bool rxn_not_mol,
      const bool molecules_not_species,
      const bool single_term
  );
  std::string generate_count_terms_for_expression(std::ostream& out, const std::string& mdl_string);

private:
  void generate_single_parameter(std::ostream& out, Json::Value& parameter);

  std::string generate_component_type(
      std::ostream& out, Json::Value& bngl_component_item, const std::string& mol_type_name);

  std::string generate_single_species_or_mol_type(
      std::ostream& out, Json::Value& molecule_list_item,
      const bool generate_species, const std::vector<std::string>& component_names = std::vector<std::string>());

  SpeciesOrMolType generate_single_species_or_mol_type_w_components(
      std::ostream& out, Json::Value& molecule_list_item);

  void generate_rxn_rule_side(std::ostream& out, Json::Value& substances_node);

  void get_surface_class_property_info(
      Json::Value& property,
      std::string& name, std::string& type_name,
      std::string& affected_mols, std::string& orientation, std::string& clamp_concentration);

  void generate_variable_rate(const std::string& rate_array_name, Json::Value& variable_rate_text);

  std::string generate_single_geometry_object(
      std::ostream& out, const int index, Json::Value& object);

  bool is_volume_mol_type(const std::string& mol_type_name);
  bool is_volume_species(const std::string& species_name);

  std::string generate_single_molecule_release_info_array(
      std::ostream& out,
      std::string& rel_site_name,
      Json::Value& release_site_list,
      Json::Value::ArrayIndex begin,
      Json::Value::ArrayIndex end);

  void generate_release_pattern(std::ostream& out, const std::string& name, std::string& delay_string);

  std::vector<std::string> get_species_to_visualize();


private:
  SharedGenData& data; // owned by MCell4Generator
  Json::Value& mcell;
  uint unnamed_surf_class_counter;
};

} /* namespace MCell */

#endif /* UTILS_DATA_MODEL_TO_PYMCELL_PYTHON_GENERATOR_H_ */
