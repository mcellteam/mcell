
#include "bng/rxn_container.h"
#include "bng/rxn_class.h"

#include <iostream>
#include <sstream>


using namespace std;

namespace BNG {

RxnContainer::~RxnContainer() {
  for (RxnClass* rc: rxn_classes) {
    delete rc;
  }

  for (RxnRule* rxn: rxn_rules) {
    delete rxn;
  }
}


void RxnContainer::reset_caches() {
  for (RxnClass* rc: rxn_classes) {
    delete rc;
  }
  rxn_classes.clear();

  species_processed_for_bimol_rxn_classes.clear();
  species_processed_for_unimol_rxn_classes.clear();
  unimol_rxn_class_map.clear();
  bimol_rxn_class_map.clear();

  // we must also erase all references to rxns classes from rxn rules
  for (RxnRule* rxn: rxn_rules) {
    rxn->reset_rxn_classes_where_used();
  }
}


void RxnContainer::update_all_mols_flags() {

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
          sp.set_flag(BNG::SPECIES_FLAG_CAN_VOLWALL);
          if (sp.id == all_mols_id || sp.id == all_vol_mols_id) {
            all_vol_mols_can_react_with_surface = true;
          }
        }
      }

      if ((sp.is_surf() || sp.id == all_mols_id)) {
        if (sp2.is_reactive_surface()) {
          if (sp.id == all_surf_mols_id || sp.id == all_mols_id) {
            all_surf_mols_can_react_with_surface = true;
          }
        }
      }
    }
  }
}


RxnClass* RxnContainer::get_or_create_empty_unimol_rxn_class(const species_id_t id) {

  auto it = unimol_rxn_class_map.find(id);

  if (it != unimol_rxn_class_map.end()) {
    return it->second;
  }
  else {
    RxnClass* new_rxn_class = new RxnClass(*this, all_species, bng_config, id);
    rxn_classes.insert(new_rxn_class);
    unimol_rxn_class_map[id] = new_rxn_class;
    return new_rxn_class;
  }
}


void RxnContainer::create_unimol_rxn_classes_for_new_species(const species_id_t new_id) {

  // find all reactions for species id
  small_vector<RxnRule*> rxns_for_new_species;
  for (RxnRule* r: rxn_rules) {
    if (r->is_unimol() && r->species_can_be_reactant(new_id, all_species)) {
      rxns_for_new_species.push_back(r);
    }
  }

  if (!rxns_for_new_species.empty()) {
    // 1) first we need to get to the instance of the reaction class for new_id
    RxnClass* rxn_class = get_or_create_empty_unimol_rxn_class(new_id);

    // create reactions classes specific for our species
    for (RxnRule* matching_rxn: rxns_for_new_species) {
      // 2) add the matching_rxn to our rxn class
      //    this also automatically updates the reaction class
      rxn_class->add_rxn_rule_no_update(matching_rxn);
    }
    rxn_class->update_rxn_pathways();

    if (bng_config.bng_verbosity_level >= 2) {
      cout << "BNG: Created a new unimolecular reaction class:\n";
      rxn_class->dump();
    }
  }

  species_processed_for_unimol_rxn_classes.insert(new_id);
}


// only for internal use
RxnClass* RxnContainer::get_or_create_empty_bimol_rxn_class(const species_id_t id1, const species_id_t id2) {

  auto it_map1 = bimol_rxn_class_map.find(id1);
  assert(it_map1 != bimol_rxn_class_map.end());

  auto it_map2 = bimol_rxn_class_map.find(id2);
  if (it_map2 == bimol_rxn_class_map.end()) {
    auto it_map2_pair = bimol_rxn_class_map.insert( make_pair(id2, SpeciesRxnClassesMap()) );
    it_map2 = it_map2_pair.first;
  }

  auto it_map1_species2 = it_map1->second.find(id2);

  // do we have a reaction class for (id1, id2)?
  if (it_map1_species2 != it_map1->second.end()) {
#ifndef NDEBUG
    auto it_map2_species1 = it_map2->second.find(id1);
    assert(it_map2_species1 != it_map2->second.end() && "Map for id1->id2 implies that there must be a map for id2->id1");
    assert(it_map1_species2->second == it_map2_species1->second && "Map for id1->id2 must be the same as a map for id2->id1");
#endif
    return it_map1_species2->second;
  }
  else {
    // create a new one
    RxnClass* new_rxn_class = new RxnClass(*this, all_species, bng_config, id1, id2);
    rxn_classes.insert(new_rxn_class);

    // insert it into maps
    it_map1->second[id2] = new_rxn_class;
    it_map2->second[id1] = new_rxn_class;
    return new_rxn_class;
  }
}


// - puts pointers to all corresponding classes to the res_classes_map
// - for bimol rxns, does not reuse already defined rxn class, e.g. when A + B was already created,
//   rxn class for B + A will be created (NOTE: might improve if needed but so far the only issue
//   are reports and printouts
// - called only from get_bimol_rxns_for_reactant
void RxnContainer::create_bimol_rxn_classes_for_new_species(const species_id_t new_id, const bool for_all_known_species) {

  // find all reactions for species id
  small_vector<RxnRule*> rxns_for_new_species;
  for (RxnRule* r: rxn_rules) {
    if (r->is_bimol() && r->species_can_be_reactant(new_id, all_species)) {
      rxns_for_new_species.push_back(r);
    }
  }

  // do nothing if this species cannot react
  if (!rxns_for_new_species.empty()) {

    // create or get rxn class map for id
    auto it = bimol_rxn_class_map.find(new_id);
    if (it == bimol_rxn_class_map.end()) {
      auto it_pair = bimol_rxn_class_map.insert( make_pair(new_id, SpeciesRxnClassesMap()) );
      it = it_pair.first;
    }

    // create reactions classes specific for our species
    // rxn_class->update_rxn_pathways may create new species, therefore we
    // must always read the current all_species contents,
    // on the other hand, polymerizing reactions might cause infinite looping therefore
    // we will limit ourselves to the species that currently exist,
    // the rxn class will be updated once a molecule of the new species (not handled here) will be created
    size_t num_species = all_species.get_species_vector().size();
    for (size_t i = 0; i < num_species; i++) {
      const Species& species = all_species.get(i);
      species_id_t second_id = species.id;

      // TODO: simplify condition - is_species_superclass check does not have to be there
      if (!for_all_known_species &&
          !species.was_instantiated()  &&
          !(new_id == second_id) && // we would miss A+A type reactions
          !all_species.is_species_superclass(second_id)) {
        // we do not care about molecules that do not exist yet, however we must process superclasses
        continue;
      }

      // don't we have a rxn class for this pair of special already?
      // (only created in different order, e.g. in A + B and now we are checking for B + A)
      auto it_class_map = bimol_rxn_class_map.find(second_id);
      if (it_class_map != bimol_rxn_class_map.end()) {
        auto it_rxn_class = it_class_map->second.find(new_id);
        if (it_rxn_class != it_class_map->second.end()) {
          // nothing to do here, this rxn class is already correctly mapped in rxn class map
          continue;
        }
      }

      small_vector<RxnRule*> applicable_rxns;

      for (RxnRule* matching_rxn: rxns_for_new_species) {

        // usually the species must be different but reactions of type A + A are allowed
        if (new_id == second_id && !matching_rxn->species_is_both_bimol_reactants(new_id, all_species)) {
          continue;
        }

        if (matching_rxn->species_can_be_bimol_reactants(new_id, second_id, all_species)) {
          // ok, we have a reaction applicable both to new_id and second_id
          // we need to add this rxn to a rxn class for these reactants

          applicable_rxns.push_back(matching_rxn);
        }
      }

      if (!applicable_rxns.empty()) {
        // get to the instance of the reaction class for (new_id, second_id) or (second_id, new_id)
        RxnClass* rxn_class = get_or_create_empty_bimol_rxn_class(new_id, second_id);

        for (RxnRule* rxn: applicable_rxns) {
          rxn_class->add_rxn_rule_no_update(rxn);
        }
        rxn_class->update_rxn_pathways();

        if (bng_config.bng_verbosity_level >= 2) {
          cout <<
              "BNG: Created or updated a new bimolecular reaction class for species " <<
              all_species.get(new_id).name << " (" << new_id << "):\n";
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

  delete rxn_class;
}


void RxnContainer::remove_unimol_rxn_class(const species_id_t id) {
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


void RxnContainer::remove_bimol_rxn_classes(const species_id_t id) {
  if (species_processed_for_bimol_rxn_classes.count(id) != 0) {
    // forget that we processed this species
    species_processed_for_bimol_rxn_classes.erase(id);

    // remove the rxn classes and their mappings for second
    // reactants

    auto it_class_map = bimol_rxn_class_map.find(id);
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
        if (reac == id) {
          continue;
        }
        auto it_reac_class_map = bimol_rxn_class_map.find(reac);
        assert(it_reac_class_map != bimol_rxn_class_map.end() &&
            "Reactant must be present in the rxn classes map");

        auto it_reac_rxn_class = it_reac_class_map->second.find(id);
        assert(it_reac_rxn_class != it_reac_class_map->second.end() && "Reverse mapping must exist");
        it_reac_class_map->second.erase(it_reac_rxn_class);
      }
    }
  }
  else {
    assert(bimol_rxn_class_map.count(id) == 0 &&
        "There must be no rxn class maps for unprocessed species");
  }
}


bool RxnContainer::has_bimol_vol_rxns() const {
  for (const BNG::RxnRule* r: rxn_rules) {
    if (r->is_bimol_vol_rxn()) {
      return true;
    }
  }
  return false;
}


void RxnContainer::dump(const bool including_rxn_rules) const {

  for (auto it_reac2: unimol_rxn_class_map) {
    cout <<
        "RxnClass for " <<
        all_species.get(it_reac2.first).name << " (" << it_reac2.first << "):\n";

    const RxnClass* rxn_class = it_reac2.second;
    assert(rxn_class != nullptr);
    rxn_class->dump("  ");
    cout << "\n";
  }

  for (auto it_reac1: bimol_rxn_class_map) {
    for (auto it_reac2: it_reac1.second) {
      cout <<
          "RxnClass for " <<
          all_species.get(it_reac1.first).name << " (" << it_reac1.first << ") + " <<
          all_species.get(it_reac2.first).name << " (" << it_reac2.first << "):\n";

      const RxnClass* rxn_class = it_reac2.second;
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
