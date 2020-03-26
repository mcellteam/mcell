/*
 * bng_engine.cpp
 *
 *  Created on: Mar 26, 2020
 *      Author: ahusar
 */

#include "bng_engine.h"

namespace BNG {

species_id_t BNGEngine::get_rxn_product_species_id(
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


CplxInstance BNGEngine::create_species_based_cplx_instance(const species_id_t id, const orientation_t o) const {
  const CplxInstance& ref = all_species.get(id);
  CplxInstance copy = ref;
  copy.set_orientation(o);
  return copy;
}


void BNGEngine::update_bimol_map_for_new_species(const species_id_t id) {
  assert(bimol_rxn_class_map.count(id) == 0 && "Not a new species");

  // find all reactions that use id as one of the reactants
  // TODO:

  // create rxn classes
  // TODO:

  // store rxn classes into the bimolecular_reactions_map
  // empty for now
  bimol_rxn_class_map[id] = SpeciesRxnClassesMap();
}




void BNGEngine::update_unimol_map_for_new_species(const species_id_t id) {
  assert(unimol_rxn_class_map.count(id) == 0 && "Not a new species");

  // find all reactions that use id as one of the reactants
  // TODO:

  // create rxn classes
  // TODO:

  // store rxn classes into the bimolecular_reactions_map
  static RxnClass empty_rxn_class;
  unimol_rxn_class_map[id] = &empty_rxn_class;
}

} // namespace BNG
