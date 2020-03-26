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

typedef std::map<species_id_t, RxnClass*> SpeciesRxnClassesMap;


typedef std::map<species_id_t, SpeciesRxnClassesMap> BimolRxnClassesMap;
typedef SpeciesRxnClassesMap UnimolRxnClassesMap;


/**
 * Owns information on reactions and species,
 * serves as a source of information for BNGEngine
 */
// TODO: better name, this will be much than a container,
// caching is done in BNGEngine
class RxnContainer {
public:
  RxnContainer(SpeciesContainer& all_species_, const BNGConfig& bng_config_)
    : all_molecules_species_id(SPECIES_ID_INVALID),
      all_volume_molecules_species_id(SPECIES_ID_INVALID),
      all_surface_molecules_species_id(SPECIES_ID_INVALID),
      all_species(all_species_),
      bng_config(bng_config_)
      {
  }

  ~RxnContainer();

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


  // this method is supposed to be used only during initialization
  void add_no_update(const RxnRule& r) {
    assert(r.is_finalized());
    // TODO: check that we don't have this rule already
    rxns.push_back(r);
  }


  const RxnClass* get_unimol_rxn_class(const species_id_t id) {
    auto it = unimol_rxn_class_map.find(id);

    // reaction maps get updated only when needed, it is not associated with addition of a new species
    // the assumption is that, after some simulation time elapsed, this will be fairly stable
    if (it == unimol_rxn_class_map.end()) {
      update_unimol_map_for_new_species(id);
      it = unimol_rxn_class_map.find(id);
    }
    return it->second;
  }



  // TODO: need orientation, check what we erased before
  const RxnClass* get_bimol_rxn_class(const species_id_t id1, const species_id_t id2) {
    assert(false);

    #if 0
    // for all reactions applicable to reacA and reacB
    BimolRxnClassesMap::const_iterator reactions_reacA_it
      = reactions.find(reacA.species_id);
    if (reactions_reacA_it == reactions.end()) {
      // no reactions at all for reacA
      return;
    }

    SpeciesRxnClassesMap::const_iterator reactions_reacA_and_reacB_it
      = reactions_reacA_it->second.find(reacB.species_id);
    if (reactions_reacA_and_reacB_it == reactions_reacA_it->second.end()) {
      return;
    }

    // there can be a single class for a unique pair of reactants,
    // TODO: check it when creating the maps
    const BNG::RxnClass* rxn_class = reactions_reacA_and_reacB_it->second;


    /* skip irrelevant reactions (i.e. non vol-surf reactions) */
    assert(rxn_class->reactants.size() == 2 && "We already checked that there must be 2 reactants");

    /* Check to see if orientation classes are zero/different */
    int test_wall = 0;
    orientation_t geomA = rxn_class->reactants[0].orientation;
    orientation_t geomB = rxn_class->reactants[1].orientation;
    if (geomA == ORIENTATION_NONE || geomB == ORIENTATION_NONE || (geomA + geomB) * (geomA - geomB) != 0) {
      matching_rxns.push_back(rxn_class);
    }
    else if (orientA != ORIENTATION_NONE && orientA * orientB * geomA * geomB > 0) {
      matching_rxns.push_back(rxn_class);
    }
    #endif
  }


  // simply looks up a reaction between 'a' and 'b',
  // this reaction must exist, asserts if not,
  // does not take species superclasses such as ALL_MOLECULES into account
  const RxnClass* get_specific_reaction_class(const species_id_t id1, const species_id_t id2) const {
    assert(false);
    /*
    assert(initialized);

    const auto& it_map_for_species = bimolecular_reactions_map.find(a.species_id);
    assert(it_map_for_species != bimolecular_reactions_map.end());

    const auto& it_res = it_map_for_species->second.find(b.species_id);
    assert(it_res != it_map_for_species->second.end());

    return it_res->second;*/
  }

  // returns null if there is no reaction for this species?
  // no -> when there is no entry in the map, this meanbs that reactants were not determined yet
  const BNG::SpeciesRxnClassesMap& get_bimol_rxns_for_reactant(const species_id_t id) {

    auto it = bimol_rxn_class_map.find(id);

    // reaction maps get updated only when needed, it is not associated with addition of a new species
    // the assumption is that, after some simulation time elapsed, this will be fairly stable
    if (it == bimol_rxn_class_map.end()) {
      update_bimol_map_for_new_species(id);
      it = bimol_rxn_class_map.find(id);
    }

    return it->second;
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

  void create_bimol_rxn_classes_for_new_species(const species_id_t id, SpeciesRxnClassesMap& res_classes_map);

private:
  RxnClass* get_or_create_empty_bimol_rxn_class(const species_id_t id1, const species_id_t id2);

  void update_unimol_map_for_new_species(const species_id_t id);
  const SpeciesRxnClassesMap& update_bimol_map_for_new_species(const species_id_t id);

private:

  // owns reaction classes
  // allocated in get_or_create_empty_bimol_rxn_class, deleted in destructor
  // the size of the vector will be changing, so we cannot take pointers to its elements
  std::vector<RxnClass*> rxn_classes;

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
  const BNGConfig& bng_config;

public:


  UnimolRxnClassesMap unimol_rxn_class_map;

  BimolRxnClassesMap bimol_rxn_class_map;

  // TODO: these should be private

  // TODO_PATHWAYS: there might be multiple reactions for 1 or 2 reactants (multiple pathways)
  //UnimolecularRxnClassesMap unimolecular_reactions_map; // created from reactions in init_simulation
  //BimolecularRxnClassesMap bimolecular_reactions_map; // created from reactions in init_simulation

};

} // namespace BNG

#endif // LIBS_BNG_RXN_CONTAINER_H_
