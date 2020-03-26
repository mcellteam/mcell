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

} // namespace BNG
