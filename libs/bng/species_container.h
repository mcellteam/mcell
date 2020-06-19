/*
 * species_container.h
 *
 *  Created on: Mar 24, 2020
 *      Author: ahusar
 */

#ifndef LIBS_BNG_SPECIES_CONTAINER_H_
#define LIBS_BNG_SPECIES_CONTAINER_H_

#include <vector>

#include "bng/bng_defines.h"
#include "bng/bng_data.h"
#include "bng/species.h"

namespace BNG {

// using templates instead of virtual methods? -> rather a template
// with virtual methods, this container would not be able to create new
// objects by its own
class SpeciesContainer {
public:
  SpeciesContainer(const BNGData& bng_data_)
    : bng_data(bng_data_),
      next_species_id(0),
      all_molecules_species_id(SPECIES_ID_INVALID),
      all_volume_molecules_species_id(SPECIES_ID_INVALID),
      all_surface_molecules_species_id(SPECIES_ID_INVALID) {
  }


  species_id_t find_or_add(const Species& new_species) {
    assert(new_species.is_finalized());

    species_id_t res = SPECIES_ID_INVALID;

    // check that this species does not exist already
    res = find(new_species);

    // add if not found
    if (res == SPECIES_ID_INVALID) {
      res = next_species_id;
      next_species_id++;

      #ifndef DEBUG
        // we do not want species with the same name
        for (const Species& s: species) {
          assert(s.name != new_species.name && "Adding species with identical name");
        }
      #endif

      species.push_back(new_species);
      species.back().id = res;
    }

    return res;
  }

  // searches for identical species
  // returns SPECIES_ID_INVALID if not found
  species_id_t find(const Species& species_to_find) {
    // simple equality comparison for now, some hashing will be needed
    for (const Species& s: species) {
      if (species_to_find.equal_ignore_id_and_flags(s)) {
        return s.id;
      }
    }
    return SPECIES_ID_INVALID;
  }

  species_id_t find(const CplxInstance& cplx_inst) {
    // simple equality comparison for now, some hashing will be needed
    for (const Species& s: species) {
      if (s.equal_cplx_instance_ignore_orientation_and_flags(cplx_inst)) {
        return s.id;
      }
    }
    return SPECIES_ID_INVALID;
  }


  species_id_t find_by_name(const std::string& name) {
    for (const Species& s: species) {
      if (s.name == name) {
        return s.id;
      }
    }
    return SPECIES_ID_INVALID;
  }

  // used from pymcell3_4
  Species* find_species_by_name(const char* name) {
    species_id_t id = find_by_name(name);
    if (id != SPECIES_ID_INVALID) {
      return &get(id);
    }
    else {
      return nullptr;
    }
  }

  species_id_t find_species_id(const CplxInstance& inst);

  species_id_t find_or_add_species(const CplxInstance& inst);

  Species& get(const species_id_t id) {
    assert(id < species.size());
    // TODO LATER: we will need some mapping, the species vector will need to be
    // 'defragmented' time from time because we cannot hold all possible species
    // in memory
    return species[id];
  }

  const Species& get(const species_id_t id) const {
    assert(id < species.size());
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

  SpeciesVector& get_species_vector() {
    return species;
  }

  void set_all_molecules_species_id(species_id_t id) {
    all_molecules_species_id = id;
  }
  void set_all_volume_molecules_species_id(species_id_t id) {
    all_volume_molecules_species_id = id;
  }
  void set_all_surface_molecules_species_id(species_id_t id) {
    all_surface_molecules_species_id = id;
  }

  species_id_t get_all_molecules_species_id() const {
    return all_molecules_species_id;
  }
  species_id_t get_all_volume_molecules_species_id() const {
    return all_volume_molecules_species_id;
  }
  species_id_t get_all_surface_molecules_species_id() const {
    return all_surface_molecules_species_id;
  }

  bool is_species_superclass(species_id_t id) const {
    return
        id == all_molecules_species_id ||
        id == all_volume_molecules_species_id ||
        id == all_surface_molecules_species_id;
  }

  void dump() const;

public:
  // FIXME: temporary, need some 'context' easily accessible by all classes
  const BNGData& bng_data;

private:
  species_id_t next_species_id;

  SpeciesVector species;

  // ids of species superclasses, SPECIES_ID_INVALID if not set
  // it might seem that this should belong into SpeciesInfo but this class needs this information
  species_id_t all_molecules_species_id;
  species_id_t all_volume_molecules_species_id;
  species_id_t all_surface_molecules_species_id;
};

} // namespace BNG

#endif /* LIBS_BNG_SPECIES_CONTAINER_H_ */
