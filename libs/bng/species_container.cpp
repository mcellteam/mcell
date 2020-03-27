#include <iostream>

#include "species_container.h"

using namespace std;

namespace BNG {

void SpeciesContainer::dump(const BNGData& bng_data) const {
  Species::dump_array(bng_data, species);
}

}
