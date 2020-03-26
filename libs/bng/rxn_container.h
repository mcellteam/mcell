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

#ifndef LIBS_BNG_RXN_CONTAINER_H_
#define LIBS_BNG_RXN_CONTAINER_H_

#include <map>

#include "bng_defines.h"
#include "rxn_rule.h"
#include "rxn_class.h"
#include "species_container.h"

namespace BNG {


/*
  TODO: caching

*/



/**
 * Owns information on reactions and species,
 * serves as a source of information for BNGEngine
 */
// TODO: maybe remove this class
class RxnContainer {
public:
  RxnContainer(SpeciesContainer& all_species_)
    : all_molecules_species_id(SPECIES_ID_INVALID),
      all_volume_molecules_species_id(SPECIES_ID_INVALID),
      all_surface_molecules_species_id(SPECIES_ID_INVALID),
      all_species(all_species_)
      {
  }

  // TODO: move to some config in bng engine
  void set_all_molecules_species_id(species_id_t id) {
    all_molecules_species_id = id;
  }
  void set_all_volume_molecules_species_id(species_id_t id) {
    all_volume_molecules_species_id = id;
  }
  void set_all_surface_molecules_species_id(species_id_t id) {
    all_surface_molecules_species_id = id;
  }


  void add(const RxnRule& r) {
    // TODO: check that we don't have this rule already
    rxns.push_back(r);
  }

#if 0

  // might return nullptr if there are none
  const SpeciesRxnClassesMap* get_specific_reactions_for_species(const species_id_t species_id) const {
    assert(false);
    /*
    if (species_id == SPECIES_ID_INVALID) {
      return nullptr;
    }

    const auto& it_map_for_species = bimolecular_reactions_map.find(species_id);
    if (it_map_for_species != bimolecular_reactions_map.end()) {
      return &it_map_for_species->second;
    }
    else {
      return nullptr;
    }*/
  }
#endif


  void dump() {
    //RxnClass::dump_array(reactions);
  }

private:

  // hold pointers to reactions
  std::vector<RxnClass> rxn_classes;

  // RxnContainer owns Rxn rules?
  // maybe just copy them after parsing
  std::vector<RxnRule> rxns;

public:
  // ids of species superclasses, SPECIES_ID_INVALID if not set
  // it might seem that this should belong into SpeciesInfo but this class needs this information
  species_id_t all_molecules_species_id;
  species_id_t all_volume_molecules_species_id;
  species_id_t all_surface_molecules_species_id;

  // owned by BNGEngine
  SpeciesContainer& all_species;

public:
  // TODO: these should be private

  // TODO_PATHWAYS: there might be multiple reactions for 1 or 2 reactants (multiple pathways)
  //UnimolecularRxnClassesMap unimolecular_reactions_map; // created from reactions in init_simulation
  //BimolecularRxnClassesMap bimolecular_reactions_map; // created from reactions in init_simulation

};

} // namespace BNG

#endif // LIBS_BNG_RXN_CONTAINER_H_
