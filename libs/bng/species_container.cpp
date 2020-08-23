#include <iostream>
#include <fstream>

#include "bng/species_container.h"

using namespace std;

namespace BNG {

species_id_t SpeciesContainer::add(const Species& new_species) {

#ifndef NDEBUG
  assert(find_full_match(new_species) == SPECIES_ID_INVALID && "Species must not exist");
  // we also don't want species with the same name
  for (const Species& s: species) {
    assert(s.name != new_species.name && "Adding species with identical name");
  }
#endif

  species_id_t res = next_species_id;
  next_species_id++;
  species.push_back(new_species);
  species.back().id = res;

  if (bng_config.debug_reactions) {
    std::cout << "BNG: Defined new species " << new_species.name << " with id " << res << "\n";
  }

  if (bng_config.reporting) {
    ofstream of;
    of.open(bng_config.get_species_report_file_name(), fstream::out | fstream::app);
    // not printing warning when file count not be opened
    if (of.is_open()) {
      of << res << ": " << new_species.to_str(bng_data) << "\n";
      of.close();
    }
  }

  return res;
}

/*
void SpeciesContainer::report_known_species() const {
  if (bng_config.reporting) {
    ofstream of;
    of.open(bng_config.get_species_report_file_name(), fstream::out | fstream::app);
    // not printing warning when file count not be opened
    if (of.is_open()) {
      for (const Species& s: species) {
        of << s.id << ": " << s.to_str(bng_data) << "\n";
      }
      of.close();
    }
  }
}*/

void SpeciesContainer::dump() const {
  Species::dump_array(bng_data, species);
}

}
