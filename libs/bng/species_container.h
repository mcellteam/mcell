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
#include "species.h"

namespace BNG {

// using templates instead of virtual methods? -> rather a template
// with virtual methods, this container would not be able to create new
// objects by its own
class SpeciesContainer {
public:
  SpeciesContainer()
    : next_species_id(0) {
  }

  species_id_t find_or_add(const Species& new_species) {
    assert(new_species.is_finalized());

    species_id_t res = SPECIES_ID_INVALID;

    // check that this species does not exist already
    if (new_species.species_id != SPECIES_ID_INVALID) {
      res = find(new_species);
    }

    // add if not found
    if (res == SPECIES_ID_INVALID) {
      res = next_species_id;
      next_species_id++;
      species.push_back(new_species);
      species.back().species_id = res;
    }

    return res;
  }

  // returns SPECIES_ID_INVALID if not found
  species_id_t find(const Species& species_to_find) {
    // simple equality comparison for now, some hashing will be needed
    for (const Species& s: species) {
      if (species_to_find.equal_except_for_id(s)) {
        return s.species_id;
      }
    }
    return SPECIES_ID_INVALID;
  }

  species_id_t find_simple_species_id(const CplxInstance& inst);


  const Species& get(const species_id_t id) const {
    assert(id < species.size());
    // TODO: we will need a mapping, the species vector will need to be
    //       defragmented time from time
    return species[id];
  }

  // for debugging
  bool is_valid_id(const species_id_t id) const {
    return id < species.size();
  }

  const CplxInstance& get_as_cplx_instance(const species_id_t id) const {
    return get(id);
  }

  orientation_t get_single_orientation(const species_id_t id) const;

  uint get_count() const {
    return species.size();
  }

  const SpeciesVector& get_species_vector() const {
    return species;
  }

  void dump(const BNGData& bng_data) const;

private:
  species_id_t next_species_id;

  SpeciesVector species;
};

} // namespace BNG

#endif /* LIBS_BNG_SPECIES_CONTAINER_H_ */
