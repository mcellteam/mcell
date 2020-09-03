#include <iostream>
#include <fstream>
#include <iomanip>

#include "bng/species_container.h"

using namespace std;

namespace BNG {

species_id_t SpeciesContainer::add(const Species& new_species) {

  Species species_copy = new_species;
  if (!species_copy.is_canonical()) {
    species_copy.canonicalize(bng_data);
  }

#ifndef NDEBUG
  assert(find_full_match(new_species) == SPECIES_ID_INVALID && "Species must not exist");
  // we also don't want species with the same name
  for (const Species& s: species) {
    assert(s.name != new_species.name && "Adding species with identical name");
  }
#endif

  species_id_t res = next_species_id;
  next_species_id++;

  // add to the species vector
  species.push_back(species_copy);
  species.back().id = res;

  // and also store canonical name for fast search
  assert(species_copy.is_canonical());
  canonical_species_map[species_copy.name] = res;

  if (bng_config.debug_reactions) {
    std::cout << "BNG: Defined new species " << species_copy.name << " with id " << res << "\n";
  }

  if (bng_config.reporting) {
    ofstream of;
    of.open(bng_config.get_species_report_file_name(), fstream::out | fstream::app);
    // not printing warning when file count not be opened
    if (of.is_open()) {
      of << res << ": " << species_copy.to_str(bng_data) <<
          ", D=" << std::setprecision(17) << species_copy.D << "\n";
      of.close();
    }
  }

  return res;
}


void SpeciesContainer::dump() const {
  Species::dump_array(bng_data, species);
}

}
