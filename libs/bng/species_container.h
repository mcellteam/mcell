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
  SpeciesContainer(const BNGData& bng_data_, const BNGConfig& bng_config_)
    : bng_data(bng_data_),
      bng_config(bng_config_),
      next_species_id(0),
      all_molecules_species_id(SPECIES_ID_INVALID),
      all_volume_molecules_species_id(SPECIES_ID_INVALID),
      all_surface_molecules_species_id(SPECIES_ID_INVALID),
      max_time_step(1.0) {
  }


  const BNGData& get_bng_data() const {
    return bng_data;
  }


  species_id_t find_or_add(Species& new_species, const bool removable = false) {
    assert(new_species.is_finalized());

    species_id_t res = SPECIES_ID_INVALID;

#ifdef DEBUG_CPLX_MATCHING
    std::cout << "Looking for new species:\n";
    new_species.dump(bng_data);
#endif

    // check that this species does not exist already
    if (!new_species.is_canonical()) {
      new_species.canonicalize();
    }
    assert(new_species.name != "");
    auto it = canonical_species_map.find(new_species.name);

    if (it == canonical_species_map.end()) {
      // add if not found
      res = add(new_species, removable);
      return res;
    }
    else {
      // return id if found
      return it->second;
    }
  }

  species_id_t add(const Species& new_species, const bool removable = false);

  void remove(const species_id_t id);

  // searches for identical species
  // returns SPECIES_ID_INVALID if not found
  species_id_t find(const Species& species_to_find) {
    // simple equality comparison for now, some hashing will be needed
    // TODO: use canonical_species_map
    for (const Species& s: species) {
      if (species_to_find.matches_fully_ignore_name_id_and_flags(s)) {
        return s.id;
      }
    }
    return SPECIES_ID_INVALID;
  }

  species_id_t find_full_match(const Cplx& cplx) const {
    // simple equality comparison for now, some hashing will be needed
    for (const Species& s: species) {
      if (s.cplx_matches_fully_ignore_orientation_and_flags(cplx)) {
        return s.id;
      }
    }
    return SPECIES_ID_INVALID;
  }


  species_id_t find_by_name(const std::string& name) const {
    // TODO: use canonical_species_map
    for (const Species& s: species) {
      if (s.name == name) {
        return s.id;
      }
    }
    return SPECIES_ID_INVALID;
  }

  // used from pymcell3_4
  const Species* find_species_by_name(const char* name) const {
    species_id_t id = find_by_name(name);
    if (id != SPECIES_ID_INVALID) {
      return &get(id);
    }
    else {
      return nullptr;
    }
  }

  Species& get(const species_id_t id) {
    assert(id != SPECIES_ID_INVALID);
    assert(id < species_id_to_index_mapping.size());
    species_index_t index = species_id_to_index_mapping[id];
    assert(index != SPECIES_INDEX_INVALID);
    assert(index < species.size());
    Species& res = species[index];
    assert(!res.is_defunct());
    return res;
  }

  const Species& get(const species_id_t id) const {
    assert(id != SPECIES_ID_INVALID);
    assert(id < species_id_to_index_mapping.size());
    species_index_t index = species_id_to_index_mapping[id];
    assert(index != SPECIES_INDEX_INVALID);
    assert(index < species.size());
    const Species& res = species[index];
    assert(!res.is_defunct());
    return res;
  }

  // for debugging
  bool is_valid_id(const species_id_t id) const {
    return
        id < species_id_to_index_mapping.size() &&
        species_id_to_index_mapping[id] < species.size();
  }

  const Cplx& get_as_cplx(const species_id_t id) const {
    return get(id);
  }

  //orientation_t get_single_orientation(const species_id_t id) const;

  uint get_count() const {
    return species.size();
  }

  const SpeciesVector& get_species_vector() const {
    return species;
  }

  SpeciesVector& get_species_vector() {
    return species;
  }
private:
  void initalize_superspecies(species_id_t id) {
    Species& sp = get(id);
    sp.set_was_instantiated(true);
    if (sp.get_num_instantiations() == 0) {
      // we want to keep the superspecies as instantiated
      sp.inc_num_instantiations();
    }
  }
public:
  void set_all_molecules_species_id(species_id_t id) {
    // superspecies are always considered to be instantiated
    all_molecules_species_id = id;
    initalize_superspecies(id);
  }
  void set_all_volume_molecules_species_id(species_id_t id) {
    all_volume_molecules_species_id = id;
    initalize_superspecies(id);
  }
  void set_all_surface_molecules_species_id(species_id_t id) {
    all_surface_molecules_species_id = id;
    initalize_superspecies(id);
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

  // flag_analyzer is a customizable class that provides interface
  // that the flags update method may query to set other flags unrelated to
  // the BNG engine itself
  void recompute_species_flags(RxnContainer& all_rxns, BaseCustomFlagsAnalyzer* flags_analyzer = nullptr) {
    for (Species& sp: species) {
      sp.update_rxn_and_custom_flags(*this, all_rxns, flags_analyzer);
    }
  }

  float_t get_max_time_step() const {
    return max_time_step;
  }

  // cleans-up the species vector by removing all species that are set as defunct
  void defragment();

  void dump() const;

public:
  const BNGData& bng_data;
  const BNGConfig& bng_config; // only debug flags is used now

private:
  species_id_t next_species_id;

  // contains mapping of molecule ids to indices to the molecules array
  std::vector<species_index_t> species_id_to_index_mapping;

  SpeciesVector species;
  std::map<std::string, species_id_t> canonical_species_map;

  // ids of species superclasses, SPECIES_ID_INVALID if not set
  // it might seem that this should belong into SpeciesInfo but this class needs this information
  species_id_t all_molecules_species_id;
  species_id_t all_volume_molecules_species_id;
  species_id_t all_surface_molecules_species_id;

  // maximal time step of any species contained in this species container,
  // this is be useful e.g. when looking for barriers in simulation
  float_t max_time_step;
};

} // namespace BNG

#endif /* LIBS_BNG_SPECIES_CONTAINER_H_ */
