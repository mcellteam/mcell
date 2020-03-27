
#include "rxn_container.h"

#include <iostream>
#include <sstream>


using namespace std;

namespace BNG {

RxnContainer::~RxnContainer() {
  for(RxnClass* rc: rxn_classes) {
    delete rc;
  }
}


void RxnContainer::create_unimol_rxn_classes_for_new_species(const species_id_t id) {
  assert(unimol_rxn_class_map.count(id) == 0 && "Not a new species");

  // find all reactions that use id as one of the reactants
  // TODO:

  // create rxn classes
  // TODO:
  //assert(false);

  // store rxn classes into the bimolecular_reactions_map
  static RxnClass empty_rxn_class(all_species, id);
  unimol_rxn_class_map[id] = &empty_rxn_class;
}

// only for internal use
RxnClass* RxnContainer::get_or_create_empty_bimol_rxn_class(const species_id_t id1, const species_id_t id2) {

  // top level rxn class maps must already exist??
  // TODO: think this through
  auto it_map1 = bimol_rxn_class_map.find(id1);
  assert(it_map1 != bimol_rxn_class_map.end());

  auto it_map2 = bimol_rxn_class_map.find(id2);
  if (it_map2 == bimol_rxn_class_map.end()) {
    SpeciesRxnClassesMap empty_map;
    bimol_rxn_class_map[id2] = empty_map;
    it_map2 = bimol_rxn_class_map.find(id2);
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
    rxn_classes.push_back(new RxnClass(all_species, id1, id2));
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
  for (RxnRule& r: rxns) {
    if (r.is_bimol() && r.species_can_be_reactant(new_id, all_species)) {
      rxns_for_new_species.push_back(&r);
    }
  }

  // 1) create rxn class map for id
  assert(bimol_rxn_class_map.count(new_id) == 0 && "Rxn class map for new species should not exist");
  SpeciesRxnClassesMap& new_rxn_class_map = bimol_rxn_class_map[new_id] = SpeciesRxnClassesMap();

  // create reactions classes specific for our species
  const SpeciesVector& species_vec = all_species.get_species_vector();
  for (const Species& s: species_vec) {

    for (RxnRule* matching_rxn: rxns_for_new_species) {

      species_id_t second_id = s.species_id;

      // usually the species must be different but reactions of type A + A are allowed
      if (new_id == second_id && !matching_rxn->species_is_both_bimol_reactants(new_id, all_species)) {
        continue;
      }

      if (matching_rxn->species_can_be_reactant(s.species_id, all_species)) {
        // ok, we have a reaction applicable both to new_id and second_id
        // we need to add this rxn to a rxn class for these reactants

        // 1) first we need to get to the instance of the reaction class for (new_id, second_id) or (second_id, new_id)
        RxnClass* rxn_class = get_or_create_empty_bimol_rxn_class(new_id, s.species_id);

        // 2) add the matching_rxn to our rxn class
        //    this also automatically updates the reaction class
        rxn_class->add_rxn_rule(bng_config, matching_rxn);
      }
    }
  }
}


species_id_t RxnContainer::get_rxn_product_species_id(
    const RxnRule* rxn, const uint product_index,
    const species_id_t reactant_a_species_id, const species_id_t reactant_b_species_id
) {
  // limited for now, no components allowed
  const CplxInstance& product = rxn->get_cplx_product(product_index);

  species_id_t res;

  if (product.is_simple()) {
    // simple species must be already defined (they are based on molecule types)
    res = all_species.find_simple_species_id(product);
    assert(res != SPECIES_ID_INVALID);
  }
  else {
    // do we have such species already or we must define a new set?
    // here are all the other tricks with
    assert(false && "TODO");
    res = SPECIES_ID_INVALID;
  }

  return res;
}


void RxnContainer::dump() const {
  set<RxnClass*> already_printed_classes;

  // TODO: unimol

  //set<RxnClass*> already_printed_classes;
  // bimol
  for (auto it_reac1: bimol_rxn_class_map) {
    for (auto it_reac2: it_reac1.second) {
      cout <<
          "RxnClass for " <<
          all_species.get(it_reac1.first).name << " (" << it_reac1.first << ") + " <<
          all_species.get(it_reac2.first).name << " (" << it_reac2.first << "):\n";

      const RxnClass* rxn_class = it_reac2.second;
      assert(rxn_class != nullptr);
      rxn_class->dump(bng_data, "  ");
      cout << "\n";
    }
  }

  cout.flush();
}

} // namespace BNG
