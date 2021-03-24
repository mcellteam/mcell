#include <iostream>
#include <fstream>
#include <iomanip>

#include "bng/species_container.h"
#include "bng/rxn_class.h"

using namespace std;

namespace BNG {

species_id_t SpeciesContainer::add(Species* new_species, const bool removable) {
  release_assert(new_species != nullptr);

#ifndef NDEBUG
  assert(find_full_match(*new_species) == SPECIES_ID_INVALID && "Species must not exist");
  // we also don't want species with the same name
  for (const Species* s: species) {
    assert(s->name != new_species->name && "Adding species with identical name");
  }
#endif

  species_id_t res = next_species_id;
  next_species_id++;

  // add to the id->index mapping and also to the species vector
  species_id_to_index_mapping.push_back(species.size());
  assert(species_id_to_index_mapping.size() == next_species_id);

  new_species->id = res;

  // and also store canonical name for fast search
  assert(new_species->is_canonical());
  canonical_species_map[new_species->name] = res;

  if (removable) {
    species.back()->set_is_removable();
  }

  // update maximal time step if needed
  if (max_time_step < new_species->time_step) {
    max_time_step = new_species->time_step;
  }

  if (bng_config.notifications.bng_verbosity_level >= 2) {
    std::cout << "BNG: Defined new species " << new_species->name << " with id " << res << "\n";
  }

  if (bng_config.rxn_and_species_report) {
    stringstream ss;
    ss <<
        res << ": " << new_species->name <<
        ", D=" << std::setprecision(17) << new_species->D <<
        ", flags: " << dynamic_cast<BaseSpeciesCplxMolFlag*>(new_species)->to_str() << "\n";

    append_to_report(bng_config.get_species_report_file_name(), ss.str());
  }

  // finally store our species
  species.push_back(new_species);

  return res;
}


void SpeciesContainer::remove(const species_id_t id) {
  // NOTE: does not remove this species from RxnContainer or anywhere else

  Species& s = get(id);
  assert(s.is_removable());

  // get compartment before we 'defunct' species because of debug checks
  compartment_id_t primary_compartment = get(id).get_primary_compartment_id();

  // do not erase directly just set that this species does not exist anymore,
  // will be physically removed on 'defragment'
  s.set_is_defunct();

  // also remove from name cache
  auto it_canonical = canonical_species_map.find(s.name);
  assert(it_canonical != canonical_species_map.end());
  canonical_species_map.erase(it_canonical);

  if (primary_compartment == COMPARTMENT_ID_NONE) {
    // species have no compartment, might need to remove it from
    // specific_to_none_compartment_species_cache
    auto it_compartment_species = none_to_specific_compartment_species_cache.find(id);
    if (it_compartment_species != none_to_specific_compartment_species_cache.end()) {
      for (auto& specific_compartment_id: it_compartment_species->second) {
        specific_to_none_compartment_species_cache.erase(specific_compartment_id.second);
      }
    }
    none_to_specific_compartment_species_cache.erase(id);
  }
  else {
    // species have compartment therefore have to be removed as a single item from
    // compartment_species_cache
    auto it_no_compartment_species = specific_to_none_compartment_species_cache.find(id);
    if (it_no_compartment_species != specific_to_none_compartment_species_cache.end()) {

      // present in cache, remove it
      assert(id == none_to_specific_compartment_species_cache[it_no_compartment_species->second][primary_compartment]);
      none_to_specific_compartment_species_cache[it_no_compartment_species->second].erase(primary_compartment);
      specific_to_none_compartment_species_cache.erase(id);
    }
  }
}


species_id_t SpeciesContainer::get_species_id_with_compartment(
    const species_id_t orig_species_id, const compartment_id_t compartment_id) {

  species_id_t no_compartment_species_id;

  compartment_id_t primary_compartment_id = get(orig_species_id).get_primary_compartment_id();
  if (primary_compartment_id != COMPARTMENT_ID_NONE) {
    // must find species that has no compartment
    auto it_map = specific_to_none_compartment_species_cache.find(orig_species_id);
    if (it_map != specific_to_none_compartment_species_cache.end()) {
      no_compartment_species_id = it_map->second;
    }
    else {
      // not cached
      // create new species with compartment_id
      Species s = get(orig_species_id);
      s.set_compartment_id(COMPARTMENT_ID_NONE);
      s.finalize_species(bng_config, false);

      // and add it
      no_compartment_species_id = find_or_add(s, true);

      // remember in cache
      none_to_specific_compartment_species_cache[no_compartment_species_id][primary_compartment_id] = orig_species_id;
      specific_to_none_compartment_species_cache[orig_species_id] = no_compartment_species_id;
    }
  }
  else {
    no_compartment_species_id = orig_species_id;
  }

  assert(is_specific_compartment_id(compartment_id));

  auto it_map = none_to_specific_compartment_species_cache.find(no_compartment_species_id);
  if (it_map != none_to_specific_compartment_species_cache.end()) {
    auto it_species = it_map->second.find(compartment_id);
    if (it_species != it_map->second.end()) {
      // found
      return it_species->second;
    }
  }

  // create new species with compartment_id
  Species s = get(no_compartment_species_id);
  s.set_compartment_id(compartment_id);
  s.finalize_species(bng_config, false);

  // and add it
  species_id_t res = find_or_add(s, true);

  // remember in cache
  none_to_specific_compartment_species_cache[no_compartment_species_id][compartment_id] = res;
  specific_to_none_compartment_species_cache[res] = no_compartment_species_id;

  return res;
}


void SpeciesContainer::defragment() {
  // based on code in MCell DefragmentationEvent::step
  size_t removed = 0;

  typedef SpeciesVector::iterator sit_t;
  sit_t it_end = species.end();

  // collect all defunct species for later

  // find first defunct molecule
  sit_t it_first_defunct =  find_if(species.begin(), it_end, [](const Species* s) -> bool { return s->is_defunct(); });
  sit_t it_copy_destination = it_first_defunct;


  while (it_first_defunct != it_end) {

    // then find the next one that is not defunct (might be it_end)
    sit_t it_next_funct = find_if(it_first_defunct, it_end, [](const Species* s) -> bool { return !s->is_defunct(); });

    // then again, find following defunct molecule
    sit_t it_second_defunct = find_if(it_next_funct, it_end, [](const Species* s) -> bool { return s->is_defunct(); });

    // items between it_first_defunct and it_next_funct will be deleted
    // we will set their ids in species_id_to_index_mapping as invalid
    for (sit_t it_update_mapping = it_first_defunct; it_update_mapping != it_next_funct; it_update_mapping++) {
      const Species* s = *it_update_mapping;
      assert(s->is_defunct());
      species_id_to_index_mapping[s->id] = SPECIES_ID_INVALID;
      delete s;
      *it_update_mapping = nullptr;      
    }

    // move data: from, to, into position
    std::copy(it_next_funct, it_second_defunct, it_copy_destination);

    removed += it_next_funct - it_first_defunct;

    // and also move destination pointer
    it_copy_destination += it_second_defunct - it_next_funct;

    it_first_defunct = it_second_defunct;
  }

  // remove everything after it_copy_destination
  if (removed != 0) {
    species.resize(species.size() - removed);
  }

  // update mapping
  size_t new_count = species.size();
  for (size_t i = 0; i < new_count; i++) {
    const Species* s = species[i];
    release_assert(s != nullptr);
    if (s->is_defunct()) {
      break;
    }
    // correct index because the species could have been moved
    species_id_to_index_mapping[s->id] = i;
  }
}


void SpeciesContainer::dump() const {
  Species::dump_array(species);
}

}
