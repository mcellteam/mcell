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

#include <string>
#include <vector>

namespace MCell {

// auxiliary struct used when generating species or molecule types
struct SpeciesOrMolType {
  SpeciesOrMolType(const std::string& name_, const bool is_species_)
    : name(name_), is_species(is_species_) {
  }

  SpeciesOrMolType(const SpeciesOrMolType& other)
    : name(other.name), is_species(other.is_species) {
  }

  bool operator == (const SpeciesOrMolType& other) const {
    return name == other.name && is_species == other.is_species;
  }

  std::string name;
  bool is_species;
};


class PythonGenerator {
public:
  PythonGenerator(Json::Value& mcell_)
    : mcell(mcell_) {
  }

  void generate_parameters(std::ostream& out);
  void generate_species_and_mol_types(std::ostream& out, std::vector<SpeciesOrMolType>& species_and_mt_info);
  void generate_geometry(std::ostream& out, std::vector<std::string>& geometry_objects);


private:
  void generate_single_parameter(std::ostream& out, Json::Value& parameter);

  std::string generate_component_type(
      std::ostream& out, Json::Value& bngl_component_item, const std::string& mol_type_name);

  std::string generate_single_species_or_mol_type(
      std::ostream& out, Json::Value& molecule_list_item,
      const bool generate_species, const std::vector<std::string>& component_names = std::vector<std::string>());

  SpeciesOrMolType generate_single_species_or_mol_type_w_components(
      std::ostream& out, Json::Value& molecule_list_item);

  std::string generate_single_geometry_object(
      std::ostream& out, const int index, Json::Value& object);


private:
  Json::Value& mcell;
};

} /* namespace MCell */

#endif /* UTILS_DATA_MODEL_TO_PYMCELL_PYTHON_GENERATOR_H_ */
