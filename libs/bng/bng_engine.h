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
