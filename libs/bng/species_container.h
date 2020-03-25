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
#include "cplx_species.h"

namespace BNG {

// using templates instead of virtual methods? -> rather a template
// with virtual methods, this container would not be able to create new
// objects by its own
// SpeciesT must be derived from CplxSpecies
template<class SpeciesT>
class SpeciesContainer {
public:
  SpeciesContainer()
    : next_species_id(0) {
  }

  species_id_t find_or_add(const SpeciesT& new_species) {
    assert(new_species.species_id == species.size());

    // check that this species does not exist already
    species_id_t res = find(new_species);

    // add if not found
    if (res == SPECIES_ID_INVALID) {
      res = next_species_id;
      next_species_id++;
      species.push_back(new_species);
    }

    return res;
  }

  // returns SPECIES_ID_INVALID if not found
  species_id_t find(const SpeciesT& species_to_find) {
    // simple equality comparison for now, some hashing will be needed
    for (const SpeciesT& s: species) {
      if (species_to_find == s) {
        return s.species_id;
      }
    }
    return SPECIES_ID_INVALID;
  }


  const SpeciesT& get(const species_id_t id) const {
    assert(id < species.size());
    // TODO: we will need a mapping, the species vector will need to be
    //       defragmented time from time
    return species[id];
  }

  const CplxInstance& get_as_cplx_instance(const species_id_t id) const {
    return get(id);
  }

  const CplxSpecies& get_as_cplx_species(const species_id_t id) const {
    return get(id);
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
  species_id_t next_species_id;

  std::vector<SpeciesT> species;
};

} // namespace BNG

#endif /* LIBS_BNG_SPECIES_CONTAINER_H_ */
