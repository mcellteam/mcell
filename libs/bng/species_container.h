/*
 * species_container.h
 *
 *  Created on: Mar 24, 2020
 *      Author: ahusar
 */

#ifndef LIBS_BNG_SPECIES_CONTAINER_H_
#define LIBS_BNG_SPECIES_CONTAINER_H_

#include <vector>

#include "bng_defines.h"

namespace BNG {

// using templates instead of virtual methods? -> rather a template
// with virtual methods, this container would not be able to create new
// objects by its own
template<class SpeciesT>
class SpeciesContainer {
public:
  SpeciesContainer()
    : next_species_id(0) {
  }

  species_id_t find_or_add(const SpeciesT& new_species) {
    assert(new_species.species_id == species.size());
    // TODO: we must check that this species does not exist already
    species_id_t res = next_species_id;
    species.push_back(new_species);
    return 0;
  }

  const SpeciesT& get(const species_id_t id) const {
    assert(id < species.size());
    // TODO: we will need a mapping, the species vector will need to be
    //       defragmented time from time
    return species[id];
  }

  uint get_count() const {
    return species.size();
  }

  const std::vector<SpeciesT>& get_species_vector() const {
    return species;
  }

  void dump() {
    //Species::dump_array(species);
  }

private:
  uint next_species_id;

  std::vector<SpeciesT> species;
};

} // namespace BNG

#endif /* LIBS_BNG_SPECIES_CONTAINER_H_ */
