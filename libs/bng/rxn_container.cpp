
#include "bng/rxn_container.h"

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
    BNG::SpeciesRxnClassesMap* rxn_classes = get_bimol_rxns_for_reactant(sp.id);
    if (rxn_classes == nullptr) {
      continue;
    }

    // go through all applicable reactants
    for (auto it: *rxn_classes) {
      const BNG::RxnClass* rxn_class = it.second;
      assert(rxn_class->is_bimol());

      const BNG::Species& sp2 = all_species.get(rxn_class->get_second_species_id(sp.id));

      // we can use is_vol/is_surf for ALL_VOLUME_MOLECULES and ALL_SURFACE_MOLECULES
      // but not for all_volume molecules because there is no
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
    rxn_classes.push_back(new RxnClass(all_species, bng_config, id));
    RxnClass* new_rxn_class = rxn_classes.back();
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
      rxn_class->add_rxn_rule(matching_rxn);
    }

    if (bng_config.debug_reactions) {
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
    rxn_classes.push_back(new RxnClass(all_species, bng_config, id1, id2));
    RxnClass* new_rxn_class = rxn_classes.back();

    // insert it into maps
    it_map1->second[id2] = new_rxn_class;
    it_map2->second[id1] = new_rxn_class;
    return new_rxn_class;
  }
}


// puts pointers to all corresponding classes to the res_classes_map
void RxnContainer::create_bimol_rxn_classes_for_new_species(const species_id_t new_id) {

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
    const SpeciesVector& species_vec = all_species.get_species_vector();
    for (const Species& s: species_vec) {

      small_vector<RxnRule*> applicable_rxns;

      for (RxnRule* matching_rxn: rxns_for_new_species) {

        species_id_t second_id = s.id;

        // usually the species must be different but reactions of type A + A are allowed
        if (new_id == second_id && !matching_rxn->species_is_both_bimol_reactants(new_id, all_species)) {
          continue;
        }

        if (matching_rxn->species_can_be_bimol_reactants(new_id, s.id, all_species)) {
          // ok, we have a reaction applicable both to new_id and second_id
          // we need to add this rxn to a rxn class for these reactants

          applicable_rxns.push_back(matching_rxn);
        }
      }

      if (!applicable_rxns.empty()) {
        // get to the instance of the reaction class for (new_id, second_id) or (second_id, new_id)
        RxnClass* rxn_class = get_or_create_empty_bimol_rxn_class(new_id, s.id);

        for (RxnRule* rxn: applicable_rxns) {
          rxn_class->add_rxn_rule(rxn);
        }

        if (bng_config.debug_reactions) {
          cout <<
              "BNG: Created or updated a new bimolecular reaction class for species " <<
              all_species.get(new_id).name << " (" << new_id << "):\n";
          rxn_class->dump("  ");
        }
      }
    } // for rxns_for_new_species
  }

  species_processed_for_bimol_rxn_classes.insert(new_id);
}


// returns true if the required reaction products were in cache
bool RxnContainer::get_cached_rxn_products(
    const RxnRule* rxn,
    const species_id_t reactant_a_species_id,
    const species_id_t reactant_b_species_id,
    std::vector<species_id_t>& res) {

  if (rxn->is_unimol()) {
    auto species_map_it = unimol_rxn_cached_product_species.find(rxn->id);
    if (species_map_it == unimol_rxn_cached_product_species.end()) {
      return false;
    }
    auto products_it = species_map_it->second.find(reactant_a_species_id);
    if (products_it == species_map_it->second.end()) {
      return false;
    }

    // found
    res = products_it->second;
    return true;
  }
  else {
    assert(rxn->is_bimol());

    auto species1_map_it = bimol_rxn_cached_product_species.find(rxn->id);
    if (species1_map_it == bimol_rxn_cached_product_species.end()) {
      return false;
    }
    auto species2_map_it = species1_map_it->second.find(reactant_a_species_id);
    if (species2_map_it == species1_map_it->second.end()) {
      return false;
    }
    auto products_it = species2_map_it->second.find(reactant_b_species_id);
    if (products_it == species2_map_it->second.end()) {
      return false;
    }

    // found
    res = products_it->second;
    return true;
  }
}


void RxnContainer::store_rxn_products_to_cache(
    const RxnRule* rxn,
    const species_id_t reactant_a_species_id,
    const species_id_t reactant_b_species_id,
    const std::vector<species_id_t>& res) {

  if (rxn->is_unimol()) {
    unimol_rxn_cached_product_species[rxn->id][reactant_a_species_id] = res;
  }
  else {
    assert(rxn->is_bimol());
    bimol_rxn_cached_product_species[rxn->id][reactant_a_species_id][reactant_b_species_id] = res;
  }
}


// this method belongs to the rxn container because reactions can have multiple
// species that it can be applied to, so caching of product species ids
// must be done here
// result must be computed on the fly because its computation may involve definition
// of new species and we must define species only when needed
// NOTE: should't this also with caching rather belong to a reaction class?
void RxnContainer::get_rxn_product_species_ids(
    const RxnRule* rxn,
    const species_id_t reactant_a_species_id,
    const species_id_t reactant_b_species_id,
    std::vector<species_id_t>& res
) {
  assert(
      (rxn->is_unimol() && reactant_b_species_id == SPECIES_ID_INVALID) ||
      (rxn->is_bimol() && (reactant_b_species_id != SPECIES_ID_INVALID || rxn->is_reactive_surface_rxn()))
  );

  // check if we did not compute this before
  // we must take the specific reactants into account
  bool found = get_cached_rxn_products(rxn, reactant_a_species_id, reactant_b_species_id, res);
  if (found) {
    return;
  }

  res.clear();

  if (rxn->is_simple()) {
    for (const CplxInstance& product: rxn->products) {
      // simple product is deterministic
      // simple species are defined mcell3 mode but in BNG mode they
      // may not have been created (they are based on molecule types)
      species_id_t id = all_species.find_full_match(product);
      if (id == SPECIES_ID_INVALID) {
        id = all_species.add(Species(product, bng_data, bng_config));
        assert(id != SPECIES_ID_INVALID);
      }
      res.push_back(id);
    }
  }
  else {
    // for a complex product, there might be multiple results
    // TODO LATER: we also need to maintain the IDs of the elementary molecules
    // but such a thing is not supported at all yet
    //
    // when the rule modifies an elementary molecule, but there is multiple
    // elementary molecules that match the reaction pattern such as in:
    // complex A(a,b~X).A(a,b~Y) and rule A(a) -> A(a!1).B(b!1).
    // there the rule can be applied to one of the distinct elementary molecules
    //
    // also having identical components inside of an elementary molecule
    // may lead to nondeterminism:
    // complex A(a,a~X) and rule A(a) -> A(a~Y),
    // there the rule can be applied to one of the distinct components

    // for now, we support reaction with a single outcome
    // and since we will depend on the input (once we will be maintaining IDs
    // of the elementary molecules), we do this computation on the fly
    vector<const CplxInstance*> reactants;
    // downcast, CplxInstance is sufficient for product computation
    reactants.push_back(dynamic_cast<const CplxInstance*>(&all_species.get(reactant_a_species_id)));
    if (rxn->is_bimol()) {
      reactants.push_back(dynamic_cast<const CplxInstance*>(&all_species.get(reactant_b_species_id)));
    }

    // we might get multiple matches on reactant(s), the numbers
    // of matches multiplied give us the total number of variants
    // a single random number then is used to choose a single variant
    // TODO: how to match this to NFSim?
    vector<CplxInstance> products;
    rxn->create_products_for_complex_rxn(
        reactants,
        products
    );

    // define the products as species
    for (const CplxInstance& product: products) {
      species_id_t id = all_species.find_or_add(Species(product, bng_data, bng_config));
      assert(id != SPECIES_ID_INVALID);
      res.push_back(id);
    }
  }

  // cache the result
  store_rxn_products_to_cache(rxn, reactant_a_species_id, reactant_b_species_id, res);
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
