/*
 * bng_engine.cpp
 *
 *  Created on: Mar 26, 2020
 *      Author: ahusar
 */
#include <iostream>
#include <sstream>

#include "bng/bng_engine.h"

using namespace std;

namespace BNG {

string BNGEngine::get_stats_report() const {
  stringstream res;

  uint num_active_species = 0;
  for (const Species& s: all_species.get_species_vector()) {
    if (s.was_instantiated()) {
      num_active_species++;
    }
  }



  res << "[" <<
      "active species " << num_active_species <<
      ", total species " << all_species.get_species_vector().size() <<
      ", total rxn classes " << all_rxns.get_num_rxn_classes() <<
      "]";
  return res.str();
}


Cplx BNGEngine::create_cplx_from_species(const species_id_t id, const orientation_t o) const {
  const Cplx& ref = all_species.get(id);
  Cplx copy = ref;
  copy.set_orientation(o);
  return copy;
}

} // namespace BNG
