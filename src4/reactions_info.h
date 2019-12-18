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

#ifndef SRC4_REACTIONS_INFO_H_
#define SRC4_REACTIONS_INFO_H_

#include "defines.h"
#include "species.h"
#include "molecule.h"
#include "reaction.h"

namespace MCell {

class Reaction;
#ifndef INDEXER_WA
typedef std::unordered_map<species_id_t, Reaction*> SpeciesReactionMap;
typedef std::unordered_map< species_id_t, SpeciesReactionMap > BimolecularReactionsMap;
#else
typedef std::map<species_id_t, Reaction*> SpeciesReactionMap;
typedef std::map<species_id_t, SpeciesReactionMap> BimolecularReactionsMap;
#endif
typedef SpeciesReactionMap UnimolecularReactionsMap;


class SpeciesInfo;

/**
 * Owns information on reactions and species,
 * mostly accessed as constant data.
 */
// TODO: move the trigger_bimolecular, trigger_bimolecular_orientation_from_mols,
// and trigger_intersect functions here?
class ReactionsInfo {
public:
  ReactionsInfo()
    : initialized(false),
      all_molecules_species_id(SPECIES_ID_INVALID),
      all_volume_molecules_species_id(SPECIES_ID_INVALID),
      all_surface_molecules_species_id(SPECIES_ID_INVALID) {
  }

  void set_all_molecules_species_id(species_id_t id) {
    all_molecules_species_id = id;
  }
  void set_all_volume_molecules_species_id(species_id_t id) {
    all_volume_molecules_species_id = id;
  }
  void set_all_surface_molecules_species_id(species_id_t id) {
    all_surface_molecules_species_id = id;
  }

  void init(const SpeciesInfo& all_species);

  bool is_initialized() const {
    return initialized;
  }

  void add(const Reaction& r) {
    reactions.push_back(r);
    initialized = false;
  }

  // simply looks up a reaction between 'a' and 'b',
  // this reaction must exist, asserts if not,
  // does not take species superclasses such as ALL_MOLECULES into account
  const Reaction* get_specific_reaction(const Molecule& a, const Molecule& b) const {
    assert(initialized);

    const auto& it_map_for_species = bimolecular_reactions_map.find(a.species_id);
    assert(it_map_for_species != bimolecular_reactions_map.end());

    const auto& it_res = it_map_for_species->second.find(b.species_id);
    assert(it_res != it_map_for_species->second.end());

    return it_res->second;
  }

  // might return nullptr if there are none
  const SpeciesReactionMap* get_specific_reactions_for_species(const species_id_t species_id) const {
    if (species_id == SPECIES_ID_INVALID) {
      return nullptr;
    }

    const auto& it_map_for_species = bimolecular_reactions_map.find(species_id);
    if (it_map_for_species != bimolecular_reactions_map.end()) {
      return &it_map_for_species->second;
    }
    else {
      return nullptr;
    }
  }

  // does not store nullptr into the resulting array
  // checks also species superclasses
  void get_all_reactions_for_reactant(const Molecule& reactant, small_vector<const SpeciesReactionMap*>& potential_reactions) const {
    potential_reactions.clear();

    // species-specific
    const SpeciesReactionMap* species_specific = get_specific_reactions_for_species(reactant.species_id);
    if (species_specific != nullptr) {
      potential_reactions.push_back(species_specific);
    }

    // all molecules
    const SpeciesReactionMap* all_molecules = get_specific_reactions_for_species(all_molecules_species_id);
    if (all_molecules != nullptr) {
      potential_reactions.push_back(all_molecules);
    }

    // all surface/volume molecules
    const SpeciesReactionMap* all_vol_surf =
        get_specific_reactions_for_species(
            reactant.is_vol() ? all_volume_molecules_species_id : all_surface_molecules_species_id);
    if (all_vol_surf != nullptr) {
      potential_reactions.push_back(all_vol_surf);
    }
  }

  void dump() {
    Reaction::dump_array(reactions);
  }

private:
  bool initialized;

  std::vector<Reaction> reactions;

  // ids of species superclasses, SPECIES_ID_INVALID if not set
  // it might seem that this should belong into SpeciesInfo but this class needs this information
  species_id_t all_molecules_species_id;
  species_id_t all_volume_molecules_species_id;
  species_id_t all_surface_molecules_species_id;

public:
  // TODO: these should be private

  // TODO_PATHWAYS: there might be multiple reactions for 1 or 2 reactants (multiple pathways)
  UnimolecularReactionsMap unimolecular_reactions_map; // created from reactions in init_simulation
  BimolecularReactionsMap bimolecular_reactions_map; // created from reactions in init_simulation

};

} // namespace mcell

#endif // SRC4_REACTIONS_INFO_H_
