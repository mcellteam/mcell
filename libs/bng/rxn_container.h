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
 * serves as a source of information for BNGEngine.
 *
 * Must be notified when a new species was created to update reaction maps.
 *
 * NOTE: Maybe rename, it is much more than a container
 */
class RxnContainer {
public:
  RxnContainer(SpeciesContainer& all_species_, const BNGData& bng_data_, const BNGConfig& bng_config_)
    : all_molecules_species_id(SPECIES_ID_INVALID),
      all_volume_molecules_species_id(SPECIES_ID_INVALID),
      all_surface_molecules_species_id(SPECIES_ID_INVALID),
      all_species(all_species_),
      bng_data(bng_data_),
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
    // TODO LATER: check that we don't have this rule already
    rxns.push_back(r);
  }


  const RxnClass* get_unimol_rxn_class(const species_id_t id) {
    auto it = unimol_rxn_class_map.find(id);

    // reaction maps get updated only when needed, it is not associated with addition of a new species
    // the assumption is that, after some simulation time elapsed, this will be fairly stable
    if (species_processed_for_unimol_rxn_classes.count(id) == 0) {
      create_unimol_rxn_classes_for_new_species(id);
      it = unimol_rxn_class_map.find(id);
    }

    if (it != unimol_rxn_class_map.end()) {
      assert(it->second != nullptr);
      assert(it->second->get_num_reactions() != 0);

      return it->second;
    }
    else {
      // no reactions for this species
      return nullptr;
    }
  }

  // simply looks up a reaction between 'a' and 'b',
  // this reaction must exist, asserts if not,
  // does not take species superclasses such as ALL_MOLECULES into account
  // order of species ids does not matter
  // get_bimol_rxn_class
  const RxnClass* get_bimol_rxn_class(const species_id_t id1, const species_id_t id2) {
    // species must exist
    assert(all_species.is_valid_id(id1));
    assert(all_species.is_valid_id(id2));

    const BNG::SpeciesRxnClassesMap* rxn_class_map_for_id1 = get_bimol_rxns_for_reactant(id1);
    if (rxn_class_map_for_id1 == nullptr) {
      // no reactions for this species at all
      return nullptr;
    }

    const auto it_res = rxn_class_map_for_id1->find(id2);

    if (it_res != rxn_class_map_for_id1->end()) {
      assert(it_res->second != nullptr);
      assert(it_res->second->get_num_reactions() != 0);

      return it_res->second;
    }
    else {
      // no reactions for this pair os species
      return nullptr;
    }
  }

  // returns null if there is no reaction for this species?
  // no -> when there is no entry in the map, this meanbs that reactants were not determined yet
  const BNG::SpeciesRxnClassesMap* get_bimol_rxns_for_reactant(const species_id_t id) {

    auto it = bimol_rxn_class_map.find(id);

    // reaction maps get updated only when needed, it is not associated with addition of a new species
    // the assumption is that, after some simulation time elapsed, this will be fairly stable
    // we must use a separate set to know whether this ID is a new one or not because
    // reaction A + B with species A also creates reaction class for B (although not a full one
    // since only reactions of B with A were considered (not B + C or or similar)
    if (species_processed_for_bimol_rxn_classes.count(id) == 0) {
      create_bimol_rxn_classes_for_new_species(id);
      it = bimol_rxn_class_map.find(id);
    }

    if (it != bimol_rxn_class_map.end()) {
      return &it->second;
    }
    else {
      // no reactions for this species at all
      return nullptr;
    }
  }

  species_id_t get_rxn_product_species_id(
      const RxnRule* rxn, const uint product_index,
      const species_id_t reactant_a_species_id, const species_id_t reactant_b_species_id
  );

  void dump() const;


private:
  RxnClass* get_or_create_empty_unimol_rxn_class(const species_id_t id);
  RxnClass* get_or_create_empty_bimol_rxn_class(const species_id_t id1, const species_id_t id2);

  void create_unimol_rxn_classes_for_new_species(const species_id_t id);

  void create_bimol_rxn_classes_for_new_species(const species_id_t id);

private:

  // owns reaction classes
  // allocated in get_or_create_empty_bimol_rxn_class, deleted in destructor
  // the size of the vector will be changing, so we cannot take pointers to its elements
  std::vector<RxnClass*> rxn_classes;

  // RxnContainer owns Rxn rules?
  // maybe just copy them after parsing
  // the size of the vector will be changing, so we cannot take pointers to its elements
  // FIXME: change to pointers
  std::vector<RxnRule> rxns;

  //
  uint_dense_hash_set<species_id_t> species_processed_for_bimol_rxn_classes;
  uint_dense_hash_set<species_id_t> species_processed_for_unimol_rxn_classes;

public:
  // TODO: move to species container
  // ids of species superclasses, SPECIES_ID_INVALID if not set
  // it might seem that this should belong into SpeciesInfo but this class needs this information
  species_id_t all_molecules_species_id;
  species_id_t all_volume_molecules_species_id;
  species_id_t all_surface_molecules_species_id;

private:
  UnimolRxnClassesMap unimol_rxn_class_map;

  BimolRxnClassesMap bimol_rxn_class_map;

  // owned by BNGEngine
  SpeciesContainer& all_species;
  const BNGData& bng_data;
  const BNGConfig& bng_config;
};

} // namespace BNG

#endif // LIBS_BNG_RXN_CONTAINER_H_
