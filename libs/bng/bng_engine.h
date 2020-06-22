/*
 * bng_engine.h
 *
 *  Created on: Jan 9, 2020
 *      Author: ahusar
 */

#ifndef LIBS_BNG_BNG_ENGINE_H_
#define LIBS_BNG_BNG_ENGINE_H_

#include "bng/bng_defines.h"
#include "bng/bng_data.h"
#include "bng/cplx_instance.h"
#include "bng/species_container.h"
#include "bng/rxn_container.h"
#include "bng/mol_type.h"
#include "bng/rxn_rule.h"

namespace BNG {

/**
 * This is the main BNG library class that owns all the data.
 *
 * One option being considered is that for parallel execution,
 * each partition in MCell will have its own BNG engine.
 * During simulation of a single time step, the contents will
 * change e.g. when a new species is created. So it must be
 * possible to synchronize the multiple BNG engines.
 */
class BNGEngine {
private:
  SpeciesContainer all_species;

  RxnContainer all_rxns;

  // data entered by user, reactions reference these data
  BNGData data;

  const BNGConfig& bng_config;

public:
  BNGEngine(const BNGConfig& bng_config_)
    : all_species(data),
      all_rxns(all_species, data, bng_config_), bng_config(bng_config_)
      {
  }

  // TODO CPLX: remove these methods

  // NOTE: use IDs for rxn patterns?
  // this function will be needed anyway
  // checks if species_id matches the reaction pattern
  bool matches(
      const CplxInstance& cplx_pattern,
      const species_id_t species_id
  ) {
    // NOTE: probably some caching will be needed
    const CplxInstance& cplx_inst = all_species.get_as_cplx_instance(species_id);
    return cplx_pattern.matches_pattern(cplx_inst);
  }

  bool matches_ignore_orientation(
      const CplxInstance& cplx_pattern,
      const species_id_t species_id
  ) {
    // NOTE: probably some caching will be needed
    const CplxInstance& cplx_inst = all_species.get_as_cplx_instance(species_id);
    return cplx_pattern.matches_pattern(cplx_inst, true);
  }

  species_id_t get_rxn_product_species_id(
      const RxnRule* rxn, const uint product_index,
      const species_id_t reactant_a_species_id, const species_id_t reactant_b_species_id
  );

  SpeciesContainer& get_all_species() { return all_species; }
  const SpeciesContainer& get_all_species() const { return all_species; }

  RxnContainer& get_all_rxns() { return all_rxns; }
  const RxnContainer& get_all_rxns() const { return all_rxns; }

  BNGData& get_data() { return data; }
  const BNGData& get_data() const { return data; }

  const BNGConfig& get_config() const { return bng_config; }

  CplxInstance create_cplx_instance_for_species(const species_id_t id, const orientation_t o) const;
};



} /* namespace BNG */

#endif /* LIBS_BNG_BNG_ENGINE_H_ */
