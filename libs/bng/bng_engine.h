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


// Each partition will have its own BNG engine,
// the contents might change a lot during execution and string comparison
// on using BNGL components we still have a way how to unify all the instances
// might need to be a template as well
//
// template - compilation time might get bad, because we need to include this every time
// need other helper classes or functions to hide the implementation
//
// SpeciesT must be derived from CplxSpecies
template<class SpeciesT>
class BNGEngine {

public:

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


  const BNG::RxnClass* get_unimol_rxn_class(const species_id_t id) {
    assert(false);
#if 0
    const UnimolecularRxnClassesMap& unimol_rxs = world->all_reactions.unimolecular_reactions_map;
    auto it = unimol_rxs.find(species_id);
    if (it == unimol_rxs.end()) {
      return nullptr;
    }
    else {
      return it->second;
    }
#endif
  }

  const BNG::RxnClass* get_bimol_rxn_class(const species_id_t id1, const species_id_t id2) {
    assert(false);

    #if 0
    // for all reactions applicable to reacA and reacB
    BimolecularRxnClassesMap::const_iterator reactions_reacA_it
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

  // does not store nullptr into the resulting array
  // checks also species superclasses
  void get_all_reactions_for_reactant(const species_id_t id, SpeciesRxnClassesMapPtrVector& potential_rxn_classes) const {
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

  // search whether two molecules can react is done
  //std::vector<ComplexSpecies> complex_species;

  BNGData& get_data() {
    return data;
  }

  const BNGData& get_data() const {
    return data;
  }

  // make private?
  SpeciesContainer<SpeciesT> all_species;


  RxnContainer all_reactions;

  // cache of complex species indices that can interact together
private:
  // data entered by user, reactions reference these data
  BNGData data;

};


template<class SpeciesT>
species_id_t BNGEngine<SpeciesT>::get_rxn_product_species_id(
    const RxnRule* rxn, const uint product_index,
    const species_id_t reactant_a_species_id, const species_id_t reactant_b_species_id
) {
  // limited for now, no components allowed
  const CplxInstance& product = rxn->get_cplx_product(product_index);

  assert(product.is_simple() && "TODO");

  // do we have such species already or we must define a new set?

  // TODO: !!! where is the list of species? -> SpeciesInfo...
  // BNG engine must be a template as well,
  //
  return SPECIES_ID_INVALID;
}

} /* namespace BNG */

#endif /* LIBS_BNG_BNG_ENGINE_H_ */
