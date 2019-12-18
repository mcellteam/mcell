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

#ifndef SRC4_REACTIONS_INFO_H_
#define SRC4_REACTIONS_INFO_H_

#include "defines.h"
#include "species.h"
#include "molecule.h"
#include "reaction.h"

namespace MCell {

class SpeciesInfo;

/**
 * Owns information on reactions and species,
 * mostly accessed as constant data.
 */
class ReactionsInfo {

public:
  ReactionsInfo()
    : initialized(false) {
  }

  void init(const SpeciesInfo& all_species);

private:

  // -------------- reaction utility methods --------------

public:

  void add(const Reaction& r) {
    reactions.push_back(r);
    initialized = false;
  }

  // TODO: fixme - doies not deal with all_molecules, etc.
  const Reaction* get_reaction(const Molecule& a, const Molecule& b) const {
    const auto& it_map_for_species = bimolecular_reactions_map.find(a.species_id);
    assert(it_map_for_species != bimolecular_reactions_map.end());
    const auto& it_res = it_map_for_species->second.find(b.species_id);
    assert(it_res != it_map_for_species->second.end());
    return it_res->second;
  }

  void dump() {
    Reaction::dump_array(reactions);
  }

private:
  bool initialized;

  std::vector<Reaction> reactions;

public:
  // TODO_PATHWAYS: there might be multiple reactions for 1 or 2 reactants (multiple pathways)
  UnimolecularReactionsMap unimolecular_reactions_map; // created from reactions in init_simulation
  BimolecularReactionsMap bimolecular_reactions_map; // created from reactions in init_simulation

};

} // namespace mcell

#endif // SRC4_REACTIONS_INFO_H_
