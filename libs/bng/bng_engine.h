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
// Does mostly caching, all the intelligence is in other classes
//
class BNGEngine {

public:

  BNGEngine(const BNGConfig& bng_config_)
    : all_rxns(all_species, bng_config_), bng_config(bng_config_)
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
  /*void can_react(
      species_id_t a,
      species_id_t b,
      small_vector<species_id_t>& reactions
  );

  void get_reactants(
      species_id_t b,
      small_vector<species_id_t>& reactions
  );
  */

  BNGData& get_data() {
    return data;
  }

  const BNGData& get_data() const {
    return data;
  }




  CplxInstance create_species_based_cplx_instance(const species_id_t id, const orientation_t o) const;


  // search whether two molecules can react is done
  //std::vector<ComplexSpecies> complex_species;


  // -------- new -----------

  /*
  // finds all reactions for species id,
  // creates reactions classes or updates existing,
  // puts pointers to all corresponding classes to the res_classes_map
  void create_rxn_classes_for_new_species(const species_id_t id, SpeciesRxnClassesMap& res_classes_map);
*/


  // make private?
  // - defintely, must be added through this engine
  SpeciesContainer all_species;

  // TODO: rename to all_rxns
  RxnContainer all_rxns;

  // cache of complex species indices that can interact together
private:
  // data entered by user, reactions reference these data
  BNGData data;

  const BNGConfig& bng_config;

};



} /* namespace BNG */

#endif /* LIBS_BNG_BNG_ENGINE_H_ */
