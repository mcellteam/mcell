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


class PythonGenerator {
public:
  PythonGenerator(Json::Value& mcell_, const std::string& output_files_prefix_, const bool bng_mode_, uint& unnamed_rxn_counter_)
    : mcell(mcell_),
      output_files_prefix(output_files_prefix_), bng_mode(bng_mode_),
      unnamed_surf_class_counter(0), unnamed_rxn_counter(unnamed_rxn_counter_) {
  }

  void generate_parameters(std::ostream& out);

  void generate_species_and_mol_types(std::ostream& out, std::vector<SpeciesOrMolType>& species_and_mt_info);
  void generate_surface_classes(std::ostream& out, std::vector<std::string>& sc_names);

  // the parameters file must be closed because we might append some code to it
  std::string generate_single_reaction_rule(std::ostream& out, Json::Value& reaction_list_item);
  void generate_reaction_rules(std::ostream& out, std::vector<IdLoc>& rxn_names);

  void generate_geometry(std::ostream& out, std::vector<std::string>& geometry_objects);

  void generate_release_sites(std::ostream& out, std::vector<std::string>& release_site_names);

  void generate_viz_outputs(std::ostream& out, const bool cellblender_viz, std::vector<std::string>& viz_output_names);


private:
  void generate_single_parameter(std::ostream& out, Json::Value& parameter);

  std::string generate_component_type(
      std::ostream& out, Json::Value& bngl_component_item, const std::string& mol_type_name);

  std::string generate_single_species_or_mol_type(
      std::ostream& out, Json::Value& molecule_list_item,
      const bool generate_species, const std::vector<std::string>& component_names = std::vector<std::string>());

  SpeciesOrMolType generate_single_species_or_mol_type_w_components(
      std::ostream& out, Json::Value& molecule_list_item);


  void get_surface_class_property_info(
      Json::Value& property,
      std::string& name, std::string& type_name, std::string& affected_mols, std::string& orientation);

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
  Json::Value& mcell;

  const std::string& output_files_prefix;
  bool bng_mode;

  uint unnamed_surf_class_counter;
  uint& unnamed_rxn_counter; // owned by MCell4Generator
};

} /* namespace MCell */

#endif /* UTILS_DATA_MODEL_TO_PYMCELL_PYTHON_GENERATOR_H_ */
