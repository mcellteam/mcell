#include <iostream>

#include "bng/species_container.h"

using namespace std;

namespace BNG {

species_id_t SpeciesContainer::find_species_id(const CplxInstance& inst) {

  // TODO LATER: we will need to use some hashing here, but let's keep it simple for now
  // see BNGEngine::matches

  for (const Species& s: species) {
    if (s.matches(inst, true)) {
      return s.id;
    }
  }
  return SPECIES_ID_INVALID;
}


void SpeciesContainer::dump() const {
  Species::dump_array(bng_data, species);
}

}
