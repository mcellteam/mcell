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

#ifndef UTILS_DATA_MODEL_TO_PYMCELL_GENERATOR_STRUCTS_H_
#define UTILS_DATA_MODEL_TO_PYMCELL_GENERATOR_STRUCTS_H_

#include <string>
#include "json/json.h"

namespace MCell {


// auxiliary struct used when generating species or molecule types
struct SpeciesOrMolType {
  SpeciesOrMolType(const std::string& name_, const bool is_species_ = true)
    : name(name_), is_species(is_species_) {
  }

  SpeciesOrMolType(const SpeciesOrMolType& other)
    : name(other.name), is_species(other.is_species) {
  }

  // ignores type
  bool operator == (const SpeciesOrMolType& other) const {
    return name == other.name;
  }

  std::string name;
  bool is_species; // mol type when false
};


struct IdLoc {
  IdLoc(const std::string& name_, const bool in_python_ = true)
    : name(name_), in_python(in_python_) {
  }

  // ignores type
  bool operator == (const IdLoc& other) const {
    return name == other.name;
  }

  std::string name;
  bool in_python; // BNGL when false
};



// data and configuration shared among generators
struct SharedGenData {
  void reset() {
    unnamed_rxn_counter = 0;
    all_species_and_mol_type_names.clear();
    all_reaction_rules_names.clear();
    bngl_reaction_rules_used_in_observables.clear();
    all_count_term_names.clear();
  }

  uint unnamed_rxn_counter;

  std::string output_files_prefix;
  bool bng_mode;
  bool debug_mode;

  std::vector<SpeciesOrMolType> all_species_and_mol_type_names;
  std::vector<IdLoc> all_reaction_rules_names;
  std::vector<std::string> bngl_reaction_rules_used_in_observables;
  std::vector<std::string> all_count_term_names;

  // mcell node of the loaded JSON file
  Json::Value mcell;
};


} // namespace MCell

#endif /* UTILS_DATA_MODEL_TO_PYMCELL_GENERATOR_STRUCTS_H_ */
