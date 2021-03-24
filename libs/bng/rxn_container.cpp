
#include "bng/rxn_container.h"
#include "bng/rxn_class.h"

#include <iostream>
#include <sstream>
#include <bitset>

using namespace std;

namespace BNG {

RxnContainer::~RxnContainer() {
  for (RxnClass* rc: rxn_classes) {
    delete rc;
  }

  for (RxnRule* rxn: rxn_rules) {
    delete rxn;
  }

  for (ReactantClass* rxn: reactant_classes_vector) {
    delete rxn;
  }
}


void RxnContainer::reset_caches() {
  for (RxnClass* rc: rxn_classes) {
    delete rc;
  }
  rxn_classes.clear();

  species_processed_for_bimol_rxn_classes.clear();
  species_processed_for_bimol_rxn_classes.shrink();
  species_processed_for_unimol_rxn_classes.clear();
  species_processed_for_bimol_rxn_classes.shrink();
  unimol_rxn_class_map.clear();
  bimol_rxn_class_map.clear();

  // we must also erase all references to rxns classes from rxn rules
  for (RxnRule* rxn: rxn_rules) {
    rxn->reset_rxn_classes_where_used();
  }
}


// this method is called only during conversion, used only for superspecies
void RxnContainer::update_all_mols_and_mol_type_compartments() {

  species_id_t all_mols_id = all_species.get_all_molecules_species_id();
  species_id_t all_vol_mols_id = all_species.get_all_volume_molecules_species_id();
  species_id_t all_surf_mols_id = all_species.get_all_surface_molecules_species_id();

  vector<species_id_t> superspecies {
    all_mols_id, all_vol_mols_id, all_surf_mols_id
  };

  // set species flags (we already have all the reactions)
  for (species_id_t species_id: superspecies) {
    BNG::Species& sp = all_species.get(species_id);

    // get reactions, this also creates all reaction classes for the species that we are processing
    // we are getting reactions for all known species because these flags must be initialized
    // before we have any instance of these species
    BNG::SpeciesRxnClassesMap* rxn_classes = get_bimol_rxns_for_reactant(sp.id, true);
    if (rxn_classes == nullptr) {
      continue;
    }

    // go through all applicable reactants
    for (auto it: *rxn_classes) {
      const BNG::RxnClass* rxn_class = it.second;
      assert(rxn_class->is_bimol());

      const BNG::Species& sp2 = all_species.get(rxn_class->get_second_species_id(sp.id));

      // we can use is_vol/is_surf for ALL_VOLUME_MOLECULES and ALL_SURFACE_MOLECULES
      // but not for all_volume molecules because there is no single flag that we can query
      if (sp.is_vol() || sp.id == all_mols_id) {
        if (sp2.is_reactive_surface()) {
          sp.set_flag(SPECIES_FLAG_CAN_VOLWALL);
          if (sp.id == all_mols_id || sp.id == all_vol_mols_id) {
            all_vol_mols_can_react_with_surface = true;
          }
        }
      }

      if ((sp.is_surf() || sp.id == all_mols_id)) {
        if (sp2.is_reactive_surface()) {
          sp.set_flag(SPECIES_FLAG_CAN_SURFWALL);
          if (sp.id == all_surf_mols_id || sp.id == all_mols_id) {
            all_surf_mols_can_react_with_surface = true;
          }
        }
      }
    }
  }
}


RxnClass* RxnContainer::get_or_create_empty_unimol_rxn_class(const species_id_t reac_id) {

  auto it = unimol_rxn_class_map.find(reac_id);

  if (it != unimol_rxn_class_map.end()) {
    return it->second;
  }
  else {
    RxnClass* new_rxn_class = new RxnClass(*this, all_species, bng_config, reac_id);
    rxn_classes.insert(new_rxn_class);
    unimol_rxn_class_map[reac_id] = new_rxn_class;
    return new_rxn_class;
  }
}


void RxnContainer::create_unimol_rxn_classes_for_new_species(const species_id_t species_id) {

  // find all reactions for species id
  small_vector<RxnRule*> rxns_for_new_species;
  for (RxnRule* r: rxn_rules) {
    if (r->is_unimol() && r->species_can_be_reactant(species_id, all_species)) {
      rxns_for_new_species.push_back(r);
    }
  }

  if (!rxns_for_new_species.empty()) {
    // 1) first we need to get to the instance of the reaction class for new_id
    RxnClass* rxn_class = get_or_create_empty_unimol_rxn_class(species_id);

    // create reactions classes specific for our species
    for (RxnRule* matching_rxn: rxns_for_new_species) {
      // 2) add the matching_rxn to our rxn class
      //    this also automatically updates the reaction class
      rxn_class->add_rxn_rule_no_update(matching_rxn);
    }
    // TODO: init rxn class pathways for debug

    if (bng_config.notifications.bng_verbosity_level >= 2) {
      cout << "BNG: Created a new unimolecular reaction class:\n";
      rxn_class->dump();
    }
  }

  species_processed_for_unimol_rxn_classes.insert(species_id);
}


RxnClass* RxnContainer::get_or_create_empty_bimol_rxn_class(const species_id_t reac1_id, const species_id_t reac2_id) {

  auto it_map1 = bimol_rxn_class_map.find(reac1_id);
  if (it_map1 == bimol_rxn_class_map.end()) {
    auto it_map1_pair = bimol_rxn_class_map.insert( make_pair(reac1_id, SpeciesRxnClassesMap()) );
    it_map1 = it_map1_pair.first;
  }

  auto it_map2 = bimol_rxn_class_map.find(reac2_id);
  if (it_map2 == bimol_rxn_class_map.end()) {
    auto it_map2_pair = bimol_rxn_class_map.insert( make_pair(reac2_id, SpeciesRxnClassesMap()) );
    it_map2 = it_map2_pair.first;
  }

  auto it_map1_species2 = it_map1->second.find(reac2_id);

  // do we have a reaction class for (id1, id2)?
  if (it_map1_species2 != it_map1->second.end()) {
#ifndef NDEBUG
    auto it_map2_species1 = it_map2->second.find(reac1_id);
    assert(it_map2_species1 != it_map2->second.end() && "Map for id1->id2 implies that there must be a map for id2->id1");
    assert(it_map1_species2->second == it_map2_species1->second && "Map for id1->id2 must be the same as a map for id2->id1");
#endif
    return it_map1_species2->second;
  }
  else {
    // create a new one
    RxnClass* new_rxn_class = new RxnClass(*this, all_species, bng_config, reac1_id, reac2_id);
    rxn_classes.insert(new_rxn_class);

    // insert it into maps
    it_map1->second[reac2_id] = new_rxn_class;
    it_map2->second[reac1_id] = new_rxn_class;
    return new_rxn_class;
  }
}


void RxnContainer::compute_reacting_classes(const ReactantClass& rc) {
  assert(rc.is_initialized());

  assert(reacting_classes.size() == rc.id);
  reacting_classes.push_back(ReactantClassIdSet());
  ReactantClassIdSet& current_set = reacting_classes.back();

  for (const ReactantClass* reacting_class: reactant_classes_set) {
    assert(reacting_class->is_initialized());

    // cross check
    // - rule_ids[0] -> reacting_class.rule_ids[1]
    // - rule_ids[1] -> reacting_class.rule_ids[0]
    bool can_react =
        (rc.reaction_id_bitsets[0] & reacting_class->reaction_id_bitsets[1]).any() ||
        (rc.reaction_id_bitsets[1] & reacting_class->reaction_id_bitsets[0]).any();

    if (can_react) {
      // mapping A -> B
      current_set.insert(reacting_class->id);

      // update also the mapping B -> A because this is a new reactant class
      reacting_classes[reacting_class->id].insert(rc.id);
    }
  }
}


// also computes reacting classes if this is a new reactant class
reactant_class_id_t RxnContainer::find_or_add_reactant_class(
    const ReactionIdBitsets& reactions_bitset_per_reactant, const bool target_only) {

  // NOTE: maybe search without creating the ReactantClass object
  ReactantClass* rc = new ReactantClass;
  rc->target_only = target_only;
  rc->reaction_id_bitsets = reactions_bitset_per_reactant;
  const auto it = reactant_classes_set.find(rc);

  if (it != reactant_classes_set.end()) {
    delete rc;
    return (*it)->id;
  }
  else {
    // add and compute reacting classes
    rc->id = next_reactant_class_id;
    next_reactant_class_id++;
    reactant_classes_set.insert(rc);
    reactant_classes_vector.push_back(rc);

    compute_reacting_classes(*rc);
    return rc->id;
  }
}


reactant_class_id_t RxnContainer::compute_reactant_class_for_species(const species_id_t species_id) {

  // prepare reactant class bitsets
  assert(rxn_rules.empty() || rxn_rules.back()->id == rxn_rules.size() - 1);
  ReactionIdBitsets reactions_bitset_per_reactant =
    { boost::dynamic_bitset<>(rxn_rules.size()), boost::dynamic_bitset<>(rxn_rules.size())};
  for (RxnRule* r: rxn_rules) {
    if (r->is_bimol_vol_rxn()) {
      std::vector<uint> indices;
      r->get_reactant_indices(species_id, all_species, indices);
      assert(indices.size() <= 2);

      if (indices.empty()) {
        continue;
      }

      if (indices.size() == 1) {
        // species matches one of reactants
        // set bit on position id
        assert(indices[0] <= 1);
        reactions_bitset_per_reactant[indices[0]][r->id] = 1;
      }
      else {
        // species matches both reactants
        assert(indices[0] + indices[1] == 1);
        reactions_bitset_per_reactant[0][r->id] = 1;
        reactions_bitset_per_reactant[1][r->id] = 1;
      }
    }
  }

  // find or add reactant class based on the computed bitsets, also compute reacting classes
  Species& s = all_species.get(species_id);
  reactant_class_id_t res = find_or_add_reactant_class(reactions_bitset_per_reactant, s.is_target_only());

#if 0
  cout <<
      reactions_bitset_per_reactant[0] << "|" << reactions_bitset_per_reactant[1] <<
      " reactant class id:" << res << ", species id:" << species_id << " (" << rxn_rules.size() << ")\n";
#endif

  return res;
}


// - puts pointers to all corresponding classes to the res_classes_map
// - for bimol rxns, does not reuse already defined rxn class, e.g. when A + B was already created,
//   rxn class for B + A will be created (NOTE: might improve if needed but so far the only issue
//   are reports and printouts
// - called only from get_bimol_rxns_for_reactant
void RxnContainer::create_bimol_rxn_classes_for_new_species(const species_id_t species_id1, const bool for_all_known_species) {

  // find all reactions for species id,
  // also define reactant class
  small_vector<RxnRule*> rxns_for_new_species;
  for (RxnRule* r: rxn_rules) {
    if (r->is_bimol() && r->species_can_be_reactant(species_id1, all_species)) {
      rxns_for_new_species.push_back(r);
    }
  }


  // do nothing if this species cannot react
  if (!rxns_for_new_species.empty()) {

    // create reactions classes specific for our species
    // rxn_class->update_rxn_pathways may create new species, therefore we
    // must always read the current all_species contents,
    // on the other hand, polymerizing reactions might cause infinite looping therefore
    // we will limit ourselves to the species that currently exist,
    // the rxn class will be updated once a molecule of the new species (not handled here) will be created
    species_index_t num_species = all_species.get_species_vector().size();
    for (species_index_t i = 0; i < num_species; i++) {
      // reading directly from the species array, not using id
      const Species* species = all_species.get_species_vector()[i];
      assert(species != nullptr);
      species_id_t species_id2 = species->id;

      // TODO: simplify condition - is_species_superclass check does not have to be there
      if (!for_all_known_species &&
          !species->was_instantiated()  &&
          !(species_id1 == species_id2) && // we would miss A+A type reactions
          !all_species.is_species_superclass(species_id2)) {
        // we do not care about molecules that do not exist yet, however we must process superclasses
        continue;
      }

      // don't we have a rxn class for this pair of special already?
      // (only created in different order, e.g. in A + B and now we are checking for B + A)
      auto it_class_map = bimol_rxn_class_map.find(species_id2);
      if (it_class_map != bimol_rxn_class_map.end()) {
        auto it_rxn_class = it_class_map->second.find(species_id1);
        if (it_rxn_class != it_class_map->second.end()) {
          // nothing to do here, this rxn class is already correctly mapped in rxn class map
          continue;
        }
      }

      small_vector<RxnRule*> applicable_rxns;

      bool has_rxn_where_both_match_both_patterns = false;

      for (RxnRule* matching_rxn: rxns_for_new_species) {

        // usually the species must be different but reactions of type A + A are allowed
        if (species_id1 == species_id2 && !matching_rxn->species_is_both_bimol_reactants(species_id1, all_species)) {
          continue;
        }

        uint reac1_index, reac2_index;
        bool both_match_both_patterns = false;
        bool reactants_match =
            matching_rxn->species_can_be_bimol_reactants(
                species_id1, species_id2, all_species,
                &reac1_index, &reac2_index, &both_match_both_patterns);

        if (both_match_both_patterns) {
          has_rxn_where_both_match_both_patterns = true;
        }

        if (reactants_match) {
          // ok, we have a reaction applicable both to new_id and second_id and compartment matches as well
          // we need to add this rxn to a rxn class for these reactants
          applicable_rxns.push_back(matching_rxn);
        }
      }

      if (!applicable_rxns.empty()) {

        species_id_t reacA_id = species_id1;
        species_id_t reacB_id = species_id2;

        // if both reactants match both patterns,
        // we must sort the reactants so that it does not matter whether the rxn class was first created
        // for reac1 or reac2,
        // otherwise, this might give us different result e.g. based on the frequency of rxn class cleanup
        // NOTE: some cases related to compartments might be missed
        if (species_id1 != species_id2 &&
            has_rxn_where_both_match_both_patterns) {
          // order by species name - the one lexicographically smaller will be the first one
          const Species& s1 = all_species.get(species_id1);
          const Species& s2 = all_species.get(species_id2);
          if (s1.name > s2.name) {
            reacA_id = species_id2;
            reacB_id = species_id1;
          }
        }

        // get to the instance of the reaction class for (new_id, second_id) or (second_id, new_id)
        RxnClass* rxn_class = get_or_create_empty_bimol_rxn_class(reacA_id, reacB_id);

        for (RxnRule* rxn: applicable_rxns) {
          rxn_class->add_rxn_rule_no_update(rxn);
        }
        // TODO: init rxn class pathways for debug

        if (bng_config.notifications.bng_verbosity_level >= 2) {
          cout <<
              "BNG: Created or updated a new bimolecular reaction class for species " <<
              all_species.get(species_id1).name << " (" << species_id2 << "):\n";
          rxn_class->dump("  ");
        }
      }

    } // for rxns_for_new_species
  }
}


void RxnContainer::delete_rxn_class(RxnClass* rxn_class) {
  // remove the pointer from the rxn classes container so that they are not deleted
  // when RxnContainer destructor is called
  assert(rxn_classes.count(rxn_class) != 0);
  rxn_classes.erase(rxn_class);

  // destructor also removes existing links from rxn rules to this rxn class
  delete rxn_class;
}


void RxnContainer::remove_unimol_rxn_classes(const species_id_t id) {
  if (species_processed_for_unimol_rxn_classes.count(id) != 0) {
    // forget that we processed this species
    species_processed_for_unimol_rxn_classes.erase(id);

    // and remove the rxn class if there is any
    auto it_rxn_class = unimol_rxn_class_map.find(id);
    if (it_rxn_class != unimol_rxn_class_map.end()) {

      delete_rxn_class(it_rxn_class->second);

      unimol_rxn_class_map.erase(it_rxn_class);
    }
  }
  else {
    assert(unimol_rxn_class_map.count(id) == 0 &&
        "There must be no rxn classes for unprocessed species");
  }
}


void RxnContainer::remove_bimol_rxn_classes(const species_id_t reac1_species_id) {

  // remove the rxn classes and their mappings for the second reactants
  if (bimol_rxn_class_map.count(reac1_species_id) != 0) {

    // forget that we processed this species (might not be present)
    species_processed_for_bimol_rxn_classes.erase(reac1_species_id);

    // remove the rxn classes and their mappings for second reactants
    auto it_class_map = bimol_rxn_class_map.find(reac1_species_id);
    if (it_class_map != bimol_rxn_class_map.end()) {

      std::vector<species_id_t> reacting_species;

      // delete all rxn classes for this species and remember
      // with which species it can react
      for (auto it_rxn_class: it_class_map->second) {
        reacting_species.push_back(it_rxn_class.first);
        delete_rxn_class(it_rxn_class.second);
      }
      bimol_rxn_class_map.erase(it_class_map);

      // ok, we removed all A + X mappings and also the rxn classes,
      // we must also remove mappings X + A
      for (species_id_t reac: reacting_species) {
        if (reac == reac1_species_id) {
          continue;
        }
        auto it_reac_class_map = bimol_rxn_class_map.find(reac);
        assert(it_reac_class_map != bimol_rxn_class_map.end() &&
            "Reactant must be present in the rxn classes map");

        auto it_reac_rxn_class = it_reac_class_map->second.find(reac1_species_id);
        assert(it_reac_rxn_class != it_reac_class_map->second.end() && "Reverse mapping must exist");
        it_reac_class_map->second.erase(it_reac_rxn_class);
      }
    }
  }
  else if (species_processed_for_bimol_rxn_classes.count(reac1_species_id) == 0) {
    // species might have been processed but if they don't react, they do not appear in the bimol_rxn_class_map
    // NOTE: not completely sure if the explanation is correct
    species_processed_for_bimol_rxn_classes.erase(reac1_species_id);
  }
}


void RxnContainer::remove_species_id_references(const species_id_t id) {
  // this method is currently used only from SpeciesCleanupEvent
  // we are assuming that there are rxn classes that use this species
  // were removed
  release_assert(species_processed_for_unimol_rxn_classes.count(id) == 0);
  release_assert(species_processed_for_bimol_rxn_classes.count(id) == 0);

  for (RxnRule* rxn: rxn_rules) {
    rxn->remove_species_id_references(id);
  }
}


void RxnContainer::remove_reactant_class(const reactant_class_id_t id) {
  for (ReactantClassIdSet& reactants: reacting_classes) {
    // remove if present
    reactants.erase(id);
  }
  // clear reactants, later, we might need to use better containers but the
  // number of reactant classes should't be high so let's keep it like this
  reacting_classes[id] = ReactantClassIdSet();

  // and delete the class itself
  ReactantClass* rc = reactant_classes_vector[id];
  assert(rc != nullptr);
  reactant_classes_set.erase(rc);
  delete rc;
  reactant_classes_vector[id] = nullptr;
}


bool RxnContainer::has_bimol_vol_rxns() const {
  for (const BNG::RxnRule* r: rxn_rules) {
    if (r->is_bimol_vol_rxn()) {
      return true;
    }
  }
  return false;
}


void RxnContainer::print_periodic_stats() const {
#define ITEM_SIZE(a) "RxnContainer: " << #a << " = " << a.size() << "\n"

  uint64_t applicable_total = 0;
  uint64_t not_applicable_total = 0;
  uint64_t rnx_class_users_total = 0;
  for (const RxnRule* rxn: rxn_rules) {
    applicable_total += rxn->species_applicable_as_any_reactant.size();
    applicable_total += rxn->species_applicable_as_reactant[0].size();
    applicable_total += rxn->species_applicable_as_reactant[1].size();

    not_applicable_total += rxn->species_not_applicable_as_any_reactant.size();
    not_applicable_total += rxn->species_not_applicable_as_reactant[0].size();
    not_applicable_total += rxn->species_not_applicable_as_reactant[1].size();
  }

  std::cout <<
      ITEM_SIZE(rxn_classes) <<
      ITEM_SIZE(species_processed_for_bimol_rxn_classes) <<
      ITEM_SIZE(species_processed_for_unimol_rxn_classes) <<
      ITEM_SIZE(unimol_rxn_class_map) <<
      ITEM_SIZE(bimol_rxn_class_map) <<
      ITEM_SIZE(rxn_rules) <<
      "RxnContainer: rxn_rules - total species applicable as reactant = " <<
        applicable_total << "\n" <<
      "RxnContainer: rxn_rules - total species not applicable as reactant = " <<
        not_applicable_total << "\n";

#undef ITEM_SIZE
}


void RxnContainer::dump(const bool including_rxn_rules) const {

  for (auto pair_species: unimol_rxn_class_map) {
    cout <<
        "RxnClass for " <<
        all_species.get(pair_species.first).name << " (" << pair_species.first << "):\n";

    const RxnClass* rxn_class = pair_species.second;
    assert(rxn_class != nullptr);
    rxn_class->dump("  ");
    cout << "\n";
  }

  for (auto pair_species1: bimol_rxn_class_map) {
    for (auto pair_species1_species2: pair_species1.second) {
      cout <<
          "RxnClass for " <<
          all_species.get(pair_species1.first).name <<
          " (" << pair_species1.first << ") + " <<
          all_species.get(pair_species1_species2.first).name <<
          " (" << pair_species1_species2.first << "):\n";

      const RxnClass* rxn_class = pair_species1_species2.second;
      assert(rxn_class != nullptr);
      rxn_class->dump("  ");
      cout << "\n";
    }
  }

  if (including_rxn_rules) {
    for (uint i = 0; i < rxn_rules.size(); i++) {
      const RxnRule* r = rxn_rules[i];
      cout << "RxnRule " << i << ": \n";
      r->dump(true, "  ");
      cout << "\n";
    }
  }

  cout.flush();
}

} // namespace BNG
