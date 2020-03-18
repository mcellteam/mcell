/******************************************************************************
 *
 * Copyright (C) 2019 by
 * The Salk Institute for Biological Studies and
 * Pittsburgh Supercomputing Center, Carnegie Mellon University
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

// FIXME: rename - this class won't contain just constants

#ifndef SRC4_SPECIES_INFO_H_
#define SRC4_SPECIES_INFO_H_

#include "defines.h"
#include "species.h"
#include "molecule.h"
#include "reaction.h"

namespace MCell {

/**
 * Owns information on reactions and species,
 * mostly accessed as constant data.
 */
class SpeciesInfo {
public:
  void add(const Species& new_species) {
    assert(new_species.species_id == species.size());
    species.push_back(new_species);
  }

  const Species& get(species_id_t id) const {
    assert(id < species.size());
    return species[id];
  }

  uint get_count() const {
    return species.size();
  }

  const std::vector<Species>& get_species_vector() const {
    return species;
  }

  void dump() {
    Species::dump_array(species);
  }

  void to_data_model(std::ostream& out) {
    out <<
        "\"define_molecules\": {\n" <<
        "\"molecule_list\": [\n";

    for (size_t i = 0; i < get_count(); i++) {
      Species &s = species[i];
      s.to_data_model(out);
      // only outputs a comma if there are more speices to come
      if(i != get_count()-1) {
        out << ",";
      }
      out << "\n";
    }
    out <<
        "],\n" <<
        "\"data_model_version\": " <<
        JSON_DM_VERSION <<
        "\n}\n";
  }

private:
  std::vector<Species> species;
};

} // namespace MCell

#endif // SRC4_SPECIES_INFO_H_
