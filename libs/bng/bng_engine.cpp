/*
 * bng_engine.cpp
 *
 *  Created on: Mar 26, 2020
 *      Author: ahusar
 */

#include "bng_engine.h"

namespace BNG {

CplxInstance BNGEngine::create_species_based_cplx_instance(const species_id_t id, const orientation_t o) const {
  const CplxInstance& ref = all_species.get(id);
  CplxInstance copy = ref;
  copy.set_orientation(o);
  return copy;
}


} // namespace BNG
