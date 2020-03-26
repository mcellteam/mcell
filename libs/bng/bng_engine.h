/*
 * bng_engine.h
 *
 *  Created on: Jan 9, 2020
 *      Author: ahusar
 */

#ifndef LIBS_BNG_BNG_ENGINE_H_
#define LIBS_BNG_BNG_ENGINE_H_

#include "bng_defines.h"
#include "bng_data.h"
#include "cplx_instance.h"
#include "species_container.h"
#include "rxn_container.h"
#include "mol_type.h"
#include "rxn_rule.h"

namespace BNG {


typedef std::map<species_id_t, RxnClass*> SpeciesRxnClassesMap;


typedef std::map<species_id_t, SpeciesRxnClassesMap> BimolRxnClassesMap;
typedef SpeciesRxnClassesMap UnimolRxnClassesMap;


// Each partition will have its own BNG engine,
// the contents might change a lot during execution and string comparison
// on using BNGL components we still have a way how to unify all the instances
// might need to be a template as well
//
// template - compilation time might get bad, because we need to include this every time
// need other helper classes or functions to hide the implementation
//
class BNGEngine {

public:

  BNGEngine()
    : all_rxns(all_species)
      {
  }

  // TODO: use IDs for rxn patterns
  // this function will be needed anyway
  // checks if species_id matches the reaction pattern
  bool matches(
      const CplxInstance& cplx_pattern,
      const species_id_t species_id
  ) {
    // TODO: caching
    const CplxInstance& cplx_inst = all_species.get_as_cplx_instance(species_id);
    return cplx_pattern.matches(cplx_inst);
  }


  bool matches_ignore_orientation(
      const CplxInstance& cplx_pattern,
      const species_id_t species_id
  ) {
    // TODO: caching
    const CplxInstance& cplx_inst = all_species.get_as_cplx_instance(species_id);
    return cplx_pattern.matches(cplx_inst, true);
  }


  species_id_t get_rxn_product_species_id(
      const RxnRule* rxn, const uint product_index,
      const species_id_t reactant_a_species_id, const species_id_t reactant_b_species_id
  );

  //component_index_t next_component_index;

  // checks cache - 2 types of caches can/cannot react
  // if not found, we need to decide
  void can_react(
      species_id_t a,
      species_id_t b,
      small_vector<species_id_t>& reactions
  );

  void get_reactants(
      species_id_t b,
      small_vector<species_id_t>& reactions
  );


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

private:
  void update_bimol_map_for_new_species(const species_id_t id);
  void update_unimol_map_for_new_species(const species_id_t id);

public:
  // does not store nullptr into the resulting array
  // checks also species superclasses
#if 0
  void get_bimol_rxns_for_reactant(const species_id_t id, SpeciesRxnClassesMapPtrVector& potential_rxn_classes) const {
    assert(false);
    /*
    potential_reactions.clear();

    // species-specific
    const SpeciesRxnClassesMap* species_specific = get_specific_reactions_for_species(reactant.species_id);
    if (species_specific != nullptr) {
      potential_reactions.push_back(species_specific);
    }

    // all molecules
    const SpeciesRxnClassesMap* all_molecules = get_specific_reactions_for_species(all_molecules_species_id);
    if (all_molecules != nullptr) {
      potential_reactions.push_back(all_molecules);
    }

    // all surface/volume molecules
    const SpeciesRxnClassesMap* all_vol_surf =
        get_specific_reactions_for_species(
            reactant.is_vol() ? all_volume_molecules_species_id : all_surface_molecules_species_id);
    if (all_vol_surf != nullptr) {
      potential_reactions.push_back(all_vol_surf);
    }
    */
  }
#endif

  CplxInstance create_species_based_cplx_instance(const species_id_t id, const orientation_t o) const;


  // search whether two molecules can react is done
  //std::vector<ComplexSpecies> complex_species;

  BNGData& get_data() {
    return data;
  }

  const BNGData& get_data() const {
    return data;
  }

  // make private?
  // - defintely, must be added through this engine
  SpeciesContainer all_species;

  // TODO: rename to all_rxns
  RxnContainer all_rxns;

  // cache of complex species indices that can interact together
private:
  // data entered by user, reactions reference these data
  BNGData data;

  UnimolRxnClassesMap unimol_rxn_class_map;

  BimolRxnClassesMap bimol_rxn_class_map;
};



} /* namespace BNG */

#endif /* LIBS_BNG_BNG_ENGINE_H_ */
