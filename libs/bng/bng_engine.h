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
#include "bng/species_container.h"
#include "bng/rxn_container.h"
#include "bng/elem_mol_type.h"
#include "bng/rxn_rule.h"
#include "bng/cplx.h"

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
    : all_species(data, bng_config_),
      all_rxns(all_species, data, bng_config_), bng_config(bng_config_)
      {
  }

  std::string get_stats_report() const;

  bool matches_ignore_orientation(
      const Cplx& cplx_pattern,
      const species_id_t species_id
  ) {
    const Cplx& cplx = all_species.get_as_cplx(species_id);
    return cplx.matches_pattern(cplx_pattern, true);
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

  Cplx create_cplx_from_species(const species_id_t id, const orientation_t o, const compartment_id_t compartment_id) const;


  // returns true if conversion was successful
  // only one compartment and volume reactions are supported now
  // returns empty string if everything went well,
  // nonempty string with error message
  std::string export_to_bngl(
      std::ostream& out_parameters,
      std::ostream& out_molecule_types,
      std::ostream& out_reaction_rules,
      const float_t volume_um3) const;

  void print_periodic_stats() const {
    std::cout << "BNG report: " << get_stats_report() << "\n";
    all_species.print_periodic_stats();
    all_rxns.print_periodic_stats();
  }

private:

  // when all_mol_types is false, no superclass species neither
  // surface classes are exported
  void export_molecule_types_as_bngl(
      std::ostream& out_parameters, std::ostream& out_molecule_types) const;
  // if error occurred, returns nonempty string with error message
  std::string export_reaction_rules_as_bngl(
      std::ostream& out_parameters, std::ostream& out_reaction_rules) const;
};



} /* namespace BNG */

#endif /* LIBS_BNG_BNG_ENGINE_H_ */
