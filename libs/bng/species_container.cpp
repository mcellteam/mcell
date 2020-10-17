#include <iostream>
#include <fstream>
#include <iomanip>

#include "bng/species_container.h"
#include "bng/rxn_class.h"

using namespace std;

namespace BNG {

species_id_t SpeciesContainer::add(const Species& new_species, const bool removable) {

  Species species_copy = new_species;
  if (!species_copy.is_canonical()) {
    species_copy.canonicalize();
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

  // add to the id->index mapping and also to the species vector
  species_id_to_index_mapping.push_back(species.size());
  assert(species_id_to_index_mapping.size() == next_species_id);

  species.push_back(species_copy);
  // TODO: clean this up, the usage of species.back() or species_copy is weird
  species.back().id = res;

  // and also store canonical name for fast search
  assert(species_copy.is_canonical());
  canonical_species_map[species_copy.name] = res;

  if (removable) {
    species.back().set_is_removable();
  }

  // update maximal time step if needed
  if (max_time_step < species_copy.time_step) {
    max_time_step = species_copy.time_step;
  }

  if (bng_config.bng_verbosity_level >= 2) {
    std::cout << "BNG: Defined new species " << species_copy.name << " with id " << res << "\n";
  }

  if (bng_config.rxn_and_species_report) {
    stringstream ss;
    ss <<
        res << ": " << species_copy.to_str() <<
        ", D=" << std::setprecision(17) << species_copy.D <<
        ", flags: " << dynamic_cast<BaseSpeciesCplxMolFlag*>(&species_copy)->to_str() << "\n";

    append_to_report(bng_config.get_species_report_file_name(), ss.str());
  }

  return res;
}


void SpeciesContainer::remove(const species_id_t id) {
  Species& s = get(id);
  assert(s.is_removable());

  // do not erase directly just set that this species does not exist anymore,
  // will be physically removed on 'defragment'
  s.set_is_defunct();

  // also remove from name cache
  auto it = canonical_species_map.find(s.name);
  assert(it != canonical_species_map.end());
  canonical_species_map.erase(it);
}


void SpeciesContainer::defragment() {
  // based on code in MCell DefragmentationEvent::step
  size_t removed = 0;

  typedef SpeciesVector::iterator sit_t;
  sit_t it_end = species.end();

  // find first defunct molecule
  sit_t it_first_defunct =  find_if(species.begin(), it_end, [](const Species & s) -> bool { return s.is_defunct(); });
  sit_t it_copy_destination = it_first_defunct;


  while (it_first_defunct != it_end) {

    // then find the next one that is not defunct (might be it_end)
    sit_t it_next_funct = find_if(it_first_defunct, it_end, [](const Species & m) -> bool { return !m.is_defunct(); });

    // then again, find following defunct molecule
    sit_t it_second_defunct = find_if(it_next_funct, it_end, [](const Species & m) -> bool { return m.is_defunct(); });

    // items between it_first_defunct and it_next_funct will be removed
    // we will set their ids in volume_molecules_id_to_index_mapping as invalid
    for (sit_t it_update_mapping = it_first_defunct; it_update_mapping != it_next_funct; it_update_mapping++) {
      const Species& vm = *it_update_mapping;
      assert(vm.is_defunct());
      species_id_to_index_mapping[vm.id] = SPECIES_ID_INVALID;
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
    species.resize(species.size() - removed, Species(bng_data));
  }

  // update mapping
  size_t new_count = species.size();
  for (size_t i = 0; i < new_count; i++) {
    const Species& s = species[i];
    if (s.is_defunct()) {
      break;
    }
    // correct index because the molecule could have been moved
    species_id_to_index_mapping[s.id] = i;
  }
}


bool SpeciesContainer::is_valid_reactant(const Reactant& reac) const {
  if (!is_valid_id(reac.species_id)) {
    return false;
  }

  set<compartment_id_t> applicable_compartments;
  get_applicable_compartments(reac.species_id, applicable_compartments);

  return applicable_compartments.count(reac.compartment_id) != 0;
}


void SpeciesContainer::get_applicable_compartments(
    const species_id_t species_id,
    std::set<compartment_id_t>& applicable_compartments) const {

  // make union from all molecule types used in this species
  applicable_compartments.clear();
  const Species& s = get(species_id);
  for (const MolInstance& mi: s.mol_instances) {
    const MolType& mt = bng_data.get_molecule_type(mi.mol_type_id);
    applicable_compartments.insert(
        mt.compartments_used_in_rxns.begin(),
        mt.compartments_used_in_rxns.end()
    );
  }
}


void SpeciesContainer::dump() const {
  Species::dump_array(bng_data, species);
}

}
