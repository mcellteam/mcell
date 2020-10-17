
#ifndef LIBS_BNG_RXN_CONTAINER_H_
#define LIBS_BNG_RXN_CONTAINER_H_

#include "bng/bng_defines.h"
#include "bng/rxn_rule.h"
#include "bng/rxn_class.h"
#include "bng/species_container.h"

namespace BNG {

// searching for a unimol rxn class:
// 1) search by species id of reactant
// 2) search by compartment of reactant (one of the precomputed options are 'any' and 'none' compartments)

typedef std::map<compartment_id_t, RxnClass*> CompartmentRxnClassesMap;
typedef CompartmentRxnClassesMap::iterator ReactantCompartmentIt;
typedef std::pair<const compartment_id_t, RxnClass*> ReactantCompartmentPair;

typedef std::map<species_id_t, CompartmentRxnClassesMap> ReactantRxnClassesMap;
typedef ReactantRxnClassesMap::iterator ReactantSpeciesIt;
typedef std::pair<const species_id_t, CompartmentRxnClassesMap> ReactantSpeciesPair;

typedef ReactantRxnClassesMap UnimolRxnClassesMap;

// searching for a bimol rxn class:
// 1) search by species id of reactant A
// 2) search by compartment of reactant A (one of the precomputed options are 'any' and 'none' compartments)
// 3) search by species id of reactant B
// 4) search by compartment of reactant B

typedef std::map<compartment_id_t, ReactantRxnClassesMap> BimolCompartmentRxnClassesMap;
typedef BimolCompartmentRxnClassesMap::iterator BimolCompartmentIt;
typedef std::pair<const compartment_id_t, ReactantRxnClassesMap> BimolCompartmentPair;

typedef std::map<species_id_t, BimolCompartmentRxnClassesMap> BimolRxnClassesMap;
typedef BimolRxnClassesMap::iterator BimolSpeciesIt;


typedef std::set<RxnClass*> RxnClassPtrSet;
typedef std::vector<RxnRule*> RxnRuleVector;

typedef std::set<species_id_t> SpeciesIdSet;

static inline BNG::RxnClass* get_rxn_class_for_any_compartment(
    const BNG::CompartmentRxnClassesMap& compartment_rxn_classes_map) {
  auto it = compartment_rxn_classes_map.find(BNG::COMPARTMENT_ID_ANY);
  assert(it != compartment_rxn_classes_map.end());
  return it->second;
}

/**
 * Owns information on reactions and species,
 * serves as a source of information for BNGEngine.
 *
 * Must be notified when a new species was created to update reaction maps.
 *
 * NOTE: Maybe rename, it is much more than a container
 */
class RxnContainer {
public:
  RxnContainer(SpeciesContainer& all_species_, BNGData& bng_data_, const BNGConfig& bng_config_)
    : all_vol_mols_can_react_with_surface(false),
      all_surf_mols_can_react_with_surface(false),
      all_species(all_species_),
      bng_data(bng_data_),
      bng_config(bng_config_) {
  }

  ~RxnContainer();

  // completely resets the rxn container, keeps only rxn rules
  void reset_caches();

  uint get_num_rxn_classes() const {
    return rxn_classes.size();
  }

  // must be called once all reactions were added or were updated
  void update_all_mols_and_mol_type_compartments();

  // this method is supposed to be used only during initialization
  rxn_rule_id_t add_and_finalize(const RxnRule& r) {
    // TODO LATER: check that we don't have this rule already

    // store a copy
    RxnRule* new_r = new RxnRule(r);
    new_r->id = rxn_rules.size();
    new_r->finalize();
    rxn_rules.push_back(new_r);
    return new_r->id;
  }

  // - might invalidate Species reference
  RxnClass* get_unimol_rxn_class(const Reactant& reac) {
    assert(all_species.is_valid_reactant(reac));

    ReactantSpeciesIt it_species = unimol_rxn_class_map.find(reac.species_id);

    // reaction maps get updated only when needed, it is not associated with addition of a new species
    // the assumption is that, after some simulation time elapsed, this will be fairly stable
    if (species_processed_for_unimol_rxn_classes.count(reac.species_id) == 0) {
      create_unimol_rxn_classes_for_new_species(reac.species_id);
      species_processed_for_unimol_rxn_classes.insert(reac.species_id);

      it_species = unimol_rxn_class_map.find(reac.species_id);
    }

    if (it_species != unimol_rxn_class_map.end()) {
      ReactantCompartmentIt it_species_comp = it_species->second.find(reac.compartment_id);
      assert(it_species_comp != it_species->second.end());

      RxnClass* res = it_species_comp->second;
      assert(res != nullptr && res->get_num_reactions() != 0);
      return res;
    }
    else {
      // no reactions for this species
      return nullptr;
    }
  }

  // frees up memory taken up by the species' rxn class that is no longer needed
  void remove_unimol_rxn_classes(const species_id_t id);

  // - simply looks up a reaction between 'a' and 'b',
  // - this reaction must exist, asserts if not,
  // - does not take species superclasses such as ALL_MOLECULES into account
  // - order of species ids does not matter
  // - might invalidate Species reference
  RxnClass* get_bimol_rxn_class(const Reactant& reac1, const Reactant& reac2) {
    // species must exist
    assert(all_species.is_valid_reactant(reac1));
    assert(all_species.is_valid_reactant(reac2));

    BNG::ReactantRxnClassesMap* ptr_species1_comp1 = get_bimol_rxns_for_reactant(reac1);
    if (ptr_species1_comp1 == nullptr) {
      // no reactions for this species at all
      return nullptr;
    }

    ReactantSpeciesIt it_species1_comp1_species2 = ptr_species1_comp1->find(reac2.species_id);

    if (it_species1_comp1_species2 != ptr_species1_comp1->end()) {

      ReactantCompartmentIt it_species1_comp1_species2_comp2 =
          it_species1_comp1_species2->second.find(reac2.compartment_id);

      if (it_species1_comp1_species2_comp2 == it_species1_comp1_species2->second.end()) {
        // no reactions for first species+compartment & second species+compartment
        return nullptr;
      }
      assert(it_species1_comp1_species2_comp2->second != nullptr);
      assert(it_species1_comp1_species2_comp2->second->get_num_reactions() != 0);

      RxnClass* res = it_species1_comp1_species2_comp2->second;
      assert(res != nullptr && res->get_num_reactions() != 0);
      return res;
    }
    else {
      // no reactions for first species+compartment & second species
      return nullptr;
    }
  }

  BNG::ReactantRxnClassesMap* get_bimol_rxns_for_reactant_any_compartment(const species_id_t species_id) {
    return get_bimol_rxns_for_reactant(Reactant(species_id, COMPARTMENT_ID_ANY));
  }

  // - returns null if there is no reaction for this species
  // - when there is no entry in the map, this means that reactions for this reactant
  //   were not determined yet and updates creates new rxn classes
  // - might invalidate Species reference
  // - if for_all_known_species is false, on rxn classes are created only for species that have 'instantiated' flag
  BNG::ReactantRxnClassesMap* get_bimol_rxns_for_reactant(const Reactant& reac, const bool for_all_known_species = false) {
    assert(all_species.is_valid_reactant(reac));

    BimolSpeciesIt it_species = bimol_rxn_class_map.find(reac.species_id);

    // did we already process this reactant?
    if (species_processed_for_bimol_rxn_classes.count(reac.species_id) == 0) {
      create_bimol_rxn_classes_for_new_species(reac.species_id, for_all_known_species);
      species_processed_for_bimol_rxn_classes.insert(reac.species_id);

      // try to find it again, maybe we did not create any rxn classes
      it_species = bimol_rxn_class_map.find(reac.species_id);
    }

    if (it_species != bimol_rxn_class_map.end()) {
      // get rxn classes for our given compartment
      auto it_species_comp = it_species->second.find(reac.compartment_id);
      if (it_species_comp == it_species->second.end()) {
        // no reactions for this species+compartment
        return nullptr;
      }

      return &it_species_comp->second;
    }
    else {
      // no reactions for this species at all
      return nullptr;
    }
  }

  // frees up memory taken up by the species' rxn classes that is no longer needed
  void remove_bimol_rxn_classes(const species_id_t reac1_species_id);

  // returns nullptr if reaction rule was not found
  RxnRule* find_rxn_rule_by_name(const std::string& name) {
    for (RxnRule* rxn_rule: rxn_rules) {
      if (rxn_rule->name == name) {
        return rxn_rule;
      }
    }
    return nullptr;
  }

  const RxnRule* get(const rxn_rule_id_t rxn_rule_id) const {
    assert(rxn_rule_id < rxn_rules.size());
    return rxn_rules[rxn_rule_id];
  }

  RxnRule* get(const rxn_rule_id_t rxn_rule_id) {
    assert(rxn_rule_id < rxn_rules.size());
    return rxn_rules[rxn_rule_id];
  }

  const RxnRuleVector& get_rxn_rules_vector() const {
    return rxn_rules;
  }

  const BNGData& get_bng_data() {
    return bng_data;
  }

  bool has_bimol_vol_rxns() const;

  void dump(const bool including_rxn_rules = false) const;

private:
  RxnClass* get_or_create_empty_unimol_rxn_class(const Reactant& reac);
  RxnClass* get_or_create_empty_bimol_rxn_class(const Reactant& reac1, const Reactant& reac2);

  void create_unimol_rxn_classes_for_new_species(const species_id_t species_id);
  void create_bimol_rxn_classes_for_new_species(const species_id_t species_id, const bool for_all_known_species);

  void delete_rxn_class(RxnClass* rxn_class);

private:
  // owns reaction classes
  // allocated in get_or_create_empty_bimol_rxn_class, deleted in destructor
  // the size of the vector will be changing, so we cannot take pointers to its elements
  // indexed by rxn_class_id_t
  RxnClassPtrSet rxn_classes;

  // RxnContainer owns Rxn rules,
  // RxnClasses use pointers to these objects
  // indexed by rxn_rule_id_t
  RxnRuleVector rxn_rules;

  // sets that remember which species were processed for rxn class generation
  SpeciesIdSet species_processed_for_bimol_rxn_classes;
  SpeciesIdSet species_processed_for_unimol_rxn_classes;

  // these two maps use reaction_hash_t keys with full orientation and
  // compartment specification
  UnimolRxnClassesMap unimol_rxn_class_map;
  BimolRxnClassesMap bimol_rxn_class_map;

  // this map allows to search for reaction classes with ignoring
  // orientation and compartment
  BimolRxnClassesMap bimol_rxn_class_any_orient_compartment_map;

public:
  // TODO: make private
  // set in update_all_mols_flags
  bool all_vol_mols_can_react_with_surface;
  bool all_surf_mols_can_react_with_surface;

private:
  // owned by BNGEngine
  SpeciesContainer& all_species;
  BNGData& bng_data;
  const BNGConfig& bng_config;
};

} // namespace BNG

#endif // LIBS_BNG_RXN_CONTAINER_H_
