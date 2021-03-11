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

typedef std::map<compartment_id_t, species_id_t> CompartmentSpeciesMap;
typedef std::map<species_id_t, CompartmentSpeciesMap> NoCompartmentToPrimaryCompartmentSpeciesMap;
typedef std::map<species_id_t, species_id_t> CompartmentToNoCompartmentSpeciesMap;


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

  ~SpeciesContainer() {
  	for (Species* s: species) {
      delete s;
    }
  }

  const BNGData& get_bng_data() const {
    return bng_data;
  }

private:
  void canonicalize_species_if_needed(Species* new_species) {
    assert(new_species->is_finalized());
    if (!new_species->is_canonical()) {
      new_species->canonicalize();
    }
    assert(new_species->name != "");
  }

public:
  // makes a copy of the passed species object
  species_id_t find_or_add(Species& new_species, const bool removable = false) {
    canonicalize_species_if_needed(&new_species);

    // check that this species does not exist already
    auto it = canonical_species_map.find(new_species.name);
    if (it == canonical_species_map.end()) {
      // make a copy and add if not found
      Species* new_species_copy = new Species(new_species);
      new_species_copy->reset_num_instantiations();
      species_id_t res = add(new_species_copy, removable);
      return res;
    }
    else {
      // return id if found
      return it->second;
    }
  }

  // SpeciesContainer takes ownership of the Species object
  // copying of species can be expensive so some variants rather use pointers
  species_id_t find_or_add_delete_if_exist(Species*& new_species, const bool removable = false) {
    canonicalize_species_if_needed(new_species);

    // check that this species does not exist already
    auto it = canonical_species_map.find(new_species->name);
    if (it == canonical_species_map.end()) {
      // make a copy and add if not found
      species_id_t res = add(new_species, removable);
      new_species = nullptr; // take over ownership
      return res;
    }
    else {
      // we already have this species, get rid of them
      delete new_species;
      new_species = nullptr;
      // and return id if found
      return it->second;
    }
  }

  species_id_t add(Species* new_species, const bool removable = false);

  void remove(const species_id_t id);

  species_id_t get_species_id_with_compartment(
      const species_id_t no_compartment_species_id, const compartment_id_t compartment_id);

  // searches for identical species
  // returns SPECIES_ID_INVALID if not found
  species_id_t find(const Species& species_to_find) {
    // simple equality comparison for now, some hashing will be needed
    // TODO: use canonical_species_map
    for (const Species* s: species) {
      if (species_to_find.matches_fully_ignore_name_id_and_flags(*s)) {
        return s->id;
      }
    }
    return SPECIES_ID_INVALID;
  }

  species_id_t find_full_match(const Cplx& cplx) const {
    // simple equality comparison for now, some hashing will be needed
    for (const Species* s: species) {
      if (s->cplx_matches_fully_ignore_orientation_and_flags(cplx)) {
        return s->id;
      }
    }
    return SPECIES_ID_INVALID;
  }

  species_id_t find_by_name(const std::string& name) const {
    // TODO: use canonical_species_map
    for (const Species* s: species) {
      if (s->name == name) {
        return s->id;
      }
    }
    return SPECIES_ID_INVALID;
  }

  Species& get(const species_id_t id) {
    assert(id != SPECIES_ID_INVALID);
    assert(id < species_id_to_index_mapping.size());
    species_index_t index = species_id_to_index_mapping[id];
    assert(index != SPECIES_INDEX_INVALID);
    assert(index < species.size());
    Species& res = *species[index];
    assert(!res.is_defunct());
    return res;
  }

  const Species& get(const species_id_t id) const {
    assert(id != SPECIES_ID_INVALID);
    assert(id < species_id_to_index_mapping.size());
    species_index_t index = species_id_to_index_mapping[id];
    assert(index != SPECIES_INDEX_INVALID);
    assert(index < species.size());
    const Species& res = *species[index];
    assert(!res.is_defunct());
    return res;
  }

  bool does_species_exist(const molecule_id_t id) {
    if (id == SPECIES_ID_INVALID) {
      return false;
    }
    if (id >= species_id_to_index_mapping.size()) {
      return false;
    }
    species_index_t index = species_id_to_index_mapping[id];
    if (index == SPECIES_INDEX_INVALID) {
      return false;
    }
    const Species* res = species[index];
    assert(res != nullptr);
    return !res->is_defunct();
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

  // to be used only in specific cases
  const std::vector<species_index_t>& get_species_id_to_index_mapping_vector() const {
    return species_id_to_index_mapping;
  }

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
  void recompute_species_flags(RxnContainer& all_rxns) {
    for (Species* sp: species) {
      release_assert(sp != nullptr);
      sp->update_rxn_and_custom_flags(*this, all_rxns);
    }
  }

  float_t get_max_time_step() const {
    return max_time_step;
  }

  // cleans-up the species vector by removing all species that are set as defunct
  void defragment();

  void print_periodic_stats() const {
    std::cout <<
        "SpeciesContainer: next_species_id = " << next_species_id << "\n" <<
        "SpeciesContainer: species_id_to_index_mapping.size() = " <<
        species_id_to_index_mapping.size() << "\n";
  }

  void dump() const;

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
  const BNGData& bng_data;
  const BNGConfig& bng_config; // only debug flags is used now

private:
  species_id_t next_species_id;

  // contains mapping of molecule ids to indices to the molecules array
  std::vector<species_index_t> species_id_to_index_mapping;

  SpeciesVector species;
  std::map<std::string, species_id_t> canonical_species_map;

  // caching of species without a compartment to species that use a single compartment for all
  // elementary molecules
  NoCompartmentToPrimaryCompartmentSpeciesMap none_to_specific_compartment_species_cache;
  // mapping in the opposite direction
  CompartmentToNoCompartmentSpeciesMap specific_to_none_compartment_species_cache;

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
