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
#include "mol_type.h"
#include "rxn_rule.h"

namespace BNG {


// Each partition will have its own BNG engine,
// the contents might change a lot during execution and string comparison
// on using BNGL components we still have a way how to unify all the instances
// might need to be a template as well
class BNGEngine {

public:

  // TODO: use IDs for rxn patterns
  // this function will be needed anyway
  // checks if species_id matches the reaction pattern
  bool matches(
      const CplxInstance& cplx_pattern,
      const species_id_t species_id
  );


  bool matches_ignore_orientation(
      const CplxInstance& cplx_pattern,
      const species_id_t species_id
  );


  species_id_t get_rxn_product_species_id(const RxnRule* rxn, const uint product_index);
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

  // cache of complex species indices that can interact together
private:
  BNGData data;
};



} /* namespace BNG */

#endif /* LIBS_BNG_BNG_ENGINE_H_ */
