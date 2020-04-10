#include <iostream>

#include "species_container.h"

using namespace std;

namespace BNG {

species_id_t SpeciesContainer::find_simple_species_id(const CplxInstance& inst) {

  // TODO LATER: we will need to use some hashing here, but let's keep it simple for now
  // see BNGEngine::matches

  for (const Species& s: species) {
    if (s.matches(inst, true)) {
      return s.species_id;
    }
  }
  return SPECIES_ID_INVALID;
}


void SpeciesContainer::dump(const BNGData& bng_data) const {
  Species::dump_array(bng_data, species);
}

}
