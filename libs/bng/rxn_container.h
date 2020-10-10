
#ifndef LIBS_BNG_RXN_CONTAINER_H_
#define LIBS_BNG_RXN_CONTAINER_H_

#include "bng/bng_defines.h"
#include "bng/rxn_rule.h"
#include "bng/rxn_class.h"
#include "bng/species_container.h"

namespace BNG {

typedef std::map<species_id_t, RxnClass*> SpeciesRxnClassesMap;

typedef std::map<species_id_t, SpeciesRxnClassesMap> BimolRxnClassesMap;
typedef SpeciesRxnClassesMap UnimolRxnClassesMap;

typedef std::set<RxnClass*> RxnClassPtrSet;
typedef std::vector<RxnRule*> RxnRuleVector;

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
  RxnContainer(SpeciesContainer& all_species_, const BNGData& bng_data_, const BNGConfig& bng_config_)
    : all_vol_mols_can_react_with_surface(false),
      all_surf_mols_can_react_with_surface(false),
      all_species(all_species_),
      bng_data(bng_data_),
      bng_config(bng_config_)
      {
  }

  ~RxnContainer();

  // completely resets the rxn container, keeps only rxn rules
  void reset_caches();

  uint get_num_rxn_classes() const {
    return rxn_classes.size();
  }

  // must be called once all reactions were added or were updated
  void update_all_mols_flags();

  // this method is supposed to be used only during initialization
  rxn_rule_id_t add_finalized_no_update(const RxnRule& r) {
    // TODO LATER: check that we don't have this rule already

    // store a copy
    RxnRule* new_r = new RxnRule(r);
    new_r->id = rxn_rules.size();
    new_r->finalize();
    rxn_rules.push_back(new_r);
    return new_r->id;
  }

  // - might invalidate Species reference
  RxnClass* get_unimol_rxn_class(const species_id_t id) {
    auto it = unimol_rxn_class_map.find(id);

    // reaction maps get updated only when needed, it is not associated with addition of a new species
    // the assumption is that, after some simulation time elapsed, this will be fairly stable
    if (species_processed_for_unimol_rxn_classes.count(id) == 0) {
      create_unimol_rxn_classes_for_new_species(id);
      it = unimol_rxn_class_map.find(id);
    }

    if (it != unimol_rxn_class_map.end()) {
      assert(it->second != nullptr);
      assert(it->second->get_num_reactions() != 0);

      return it->second;
    }
    else {
      // no reactions for this species
      return nullptr;
    }
  }

  // frees up memory taken up by the species' rxn class that is no longer needed
  void remove_unimol_rxn_class(const species_id_t id);

  // - simply looks up a reaction between 'a' and 'b',
  // - this reaction must exist, asserts if not,
  // - does not take species superclasses such as ALL_MOLECULES into account
  // - order of species ids does not matter
  // - might invalidate Species reference
  RxnClass* get_bimol_rxn_class(const species_id_t id1, const species_id_t id2) {
    // species must exist
    assert(all_species.is_valid_id(id1));
    assert(all_species.is_valid_id(id2));

    BNG::SpeciesRxnClassesMap* rxn_class_map_for_id1 = get_bimol_rxns_for_reactant(id1);
    if (rxn_class_map_for_id1 == nullptr) {
      // no reactions for this species at all
      return nullptr;
    }

    const auto it_res = rxn_class_map_for_id1->find(id2);

    if (it_res != rxn_class_map_for_id1->end()) {
      assert(it_res->second != nullptr);
      assert(it_res->second->get_num_reactions() != 0);

      return it_res->second;
    }
    else {
      // no reactions for this pair os species
      return nullptr;
    }
  }

  // - returns null if there is no reaction for this species
  // - when there is no entry in the map, this means that reactions for this reactant
  //   were not determined yet and updates creates new rxn classes
  // - might invalidate Species reference
  BNG::SpeciesRxnClassesMap* get_bimol_rxns_for_reactant(const species_id_t id, const bool for_all_known_species = false) {

    auto it = bimol_rxn_class_map.find(id);

    // reaction maps get updated only when needed, it is not associated with addition of a new species
    // the assumption is that, after some simulation time elapsed, this will be fairly stable
    // we must use a separate set to know whether this ID is a new one or not because
    // reaction A + B with species A also creates reaction class for B (although not a full one
    // since only reactions of B with A were considered (not B + C or or similar)
    // ??? comment may not be up to date anymore
    /// TODO add   it != bimol_rxn_class_map.end() &&
    if (species_processed_for_bimol_rxn_classes.count(id) == 0) {
      create_bimol_rxn_classes_for_new_species(id, for_all_known_species);
      species_processed_for_bimol_rxn_classes.insert(id);

      // try to find it again, maybe we did not create any rxn classes
      it = bimol_rxn_class_map.find(id);
    }

    if (it != bimol_rxn_class_map.end()) {
      return &it->second;
    }
    else {
      // no reactions for this species at all
      return nullptr;
    }
  }

  // frees up memory taken up by the species' rxn classes that is no longer needed
  void remove_bimol_rxn_classes(const species_id_t id);

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
  RxnClass* get_or_create_empty_unimol_rxn_class(const species_id_t id);
  RxnClass* get_or_create_empty_bimol_rxn_class(const species_id_t id1, const species_id_t id2);

  void create_unimol_rxn_classes_for_new_species(const species_id_t id);
  void create_bimol_rxn_classes_for_new_species(const species_id_t id, const bool for_all_known_species);

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
  uint_dense_hash_set<species_id_t> species_processed_for_bimol_rxn_classes;
  uint_dense_hash_set<species_id_t> species_processed_for_unimol_rxn_classes;

  UnimolRxnClassesMap unimol_rxn_class_map;

  BimolRxnClassesMap bimol_rxn_class_map;

public:
  // TODO: make private
  // set in update_all_mols_flags
  bool all_vol_mols_can_react_with_surface;
  bool all_surf_mols_can_react_with_surface;

private:
  // owned by BNGEngine
  SpeciesContainer& all_species;
  const BNGData& bng_data;
  const BNGConfig& bng_config;
};

} // namespace BNG

#endif // LIBS_BNG_RXN_CONTAINER_H_
