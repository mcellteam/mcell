
#ifndef LIBS_BNG_RXN_CONTAINER_H_
#define LIBS_BNG_RXN_CONTAINER_H_

#define BOOST_ALLOW_DEPRECATED_HEADERS
#include "libs/boost/dynamic_bitset.hpp"

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

typedef uint_set<species_id_t> SpeciesIdSet;

// used always with two elements, the contained bitsets have size of
// rxn_rules.size(), first one is for the first reactant, second for the second reactant
// WARNING: the number of reaction rules cannot change in runtime for this to work
typedef std::vector<boost::dynamic_bitset<>> ReactionIdBitsets;


class ReactantClass {
public:
  bool is_initialized() const {
    return reaction_id_bitsets.size() == 2;
  }

  // used for lookup, does not compare id
  bool operator < (const ReactantClass& other) const {
    assert(is_initialized());
    assert(other.is_initialized());

    if (target_only != other.target_only) {
      return (target_only?1:0) < (other.target_only?1:0);
    }
    else if (reaction_id_bitsets[0] != other.reaction_id_bitsets[0]) {
      return reaction_id_bitsets[0] < other.reaction_id_bitsets[0];
    }
    else {
      return reaction_id_bitsets[1] < other.reaction_id_bitsets[1];
    }
  }

  reactant_class_id_t id;
  bool target_only;
  ReactionIdBitsets reaction_id_bitsets;
};

typedef std::vector<ReactantClass*> ReactantClassesVector;

class ReactantClassLessPtr {
public:
  bool operator() (const ReactantClass* a, const ReactantClass* b) const {
    assert(a != nullptr && b != nullptr);
    return *a < *b;
  }
};


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
    : next_reactant_class_id(0),
      all_vol_mols_can_react_with_surface(false),
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
  // - returns nullptr when there are no rxns, never returns an empty rxn class
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
  void remove_unimol_rxn_classes(const species_id_t id);

  // - simply looks up a reaction between 'a' and 'b',
  // - this reaction must exist, asserts if not,
  // - does not take species superclasses such as ALL_MOLECULES into account
  // - order of species ids does not matter
  // - might invalidate Species reference
  // - returns nullptr when there are no rxns, never returns an empty rxn class
  RxnClass* get_bimol_rxn_class(const species_id_t reac1_id, const species_id_t reac2_id) {
    // species must exist
    assert(all_species.is_valid_id(reac1_id));
    assert(all_species.is_valid_id(reac2_id));

    BNG::SpeciesRxnClassesMap* rxn_class_map_for_id1 = get_bimol_rxns_for_reactant(reac1_id);
    if (rxn_class_map_for_id1 == nullptr) {
      // no reactions for this species at all
      return nullptr;
    }

    const auto it_res = rxn_class_map_for_id1->find(reac2_id);

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
  // - if for_all_known_species is false, on rxn classes are created only for species that have 'instantiated' flag
  BNG::SpeciesRxnClassesMap* get_bimol_rxns_for_reactant(const species_id_t reac_id, const bool for_all_known_species = false) {
    assert(all_species.is_valid_id(reac_id));

    auto it = bimol_rxn_class_map.find(reac_id);

    // reaction maps get updated only when needed, it is not associated with addition of a new species
    // the assumption is that, after some simulation time elapsed, this will be fairly stable
    // we must use a separate set to know whether this ID is a new one or not because
    // reaction A + B with species A also creates reaction class for B (although not a full one
    // since only reactions of B with A were considered (not B + C or or similar)
    // ??? comment may not be up to date anymore
    /// TODO add   it != bimol_rxn_class_map.end() &&
    if (species_processed_for_bimol_rxn_classes.count(reac_id) == 0) {
      create_bimol_rxn_classes_for_new_species(reac_id, for_all_known_species);
      species_processed_for_bimol_rxn_classes.insert(reac_id);

      // try to find it again, maybe we did not create any rxn classes
      it = bimol_rxn_class_map.find(reac_id);
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
  void remove_bimol_rxn_classes(const species_id_t reac1_species_id);

  void remove_species_id_references(const species_id_t id);

  const ReactantClassIdSet& get_reacting_classes(const species_id_t species_id) {

    // get or compute species reactant class
    reactant_class_id_t reactant_class_id;
    Species& species = all_species.get(species_id);
    assert(species.is_vol() && "Reacting classes are currently supported only for volume molecules");

    if (species.has_valid_reactant_class_id()) {
      reactant_class_id = species.get_reactant_class_id();
    }
    else {
      reactant_class_id = compute_reactant_class_for_species(species_id);
      species.set_reactant_class_id(reactant_class_id);
    }

    assert(reactant_class_id < reacting_classes.size());
    return reacting_classes[reactant_class_id];
  }

  const ReactantClass& get_reactant_class(const reactant_class_id_t id) {
    assert(id != REACTANT_CLASS_ID_INVALID);
    assert(id < reactant_classes_vector.size());
    return *reactant_classes_vector[id];
  }

  size_t get_num_existing_reactant_classes() const {
#ifdef NDEBUG
    size_t num = 0;
    for (const auto& rc: reactant_classes_vector) {
      if (rc != nullptr) {
        num++;
      }
    }
    assert(num == reactant_classes_set.size());
#endif
    return reactant_classes_set.size();
  }

  // may contain nullptr items
  ReactantClassesVector get_reactant_classes() const {
    return reactant_classes_vector;
  }

  void remove_reactant_class(const reactant_class_id_t id);

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

  void print_periodic_stats() const;

  void dump(const bool including_rxn_rules = false) const;

private:
  RxnClass* get_or_create_empty_unimol_rxn_class(const species_id_t reac_id);
  RxnClass* get_or_create_empty_bimol_rxn_class(const species_id_t reac1_id, const species_id_t reac2_id);

  void create_unimol_rxn_classes_for_new_species(const species_id_t species_id);
  void create_bimol_rxn_classes_for_new_species(const species_id_t species_id, const bool for_all_known_species);

  void delete_rxn_class(RxnClass* rxn_class);

  void compute_reacting_classes(const ReactantClass& rc);
  reactant_class_id_t find_or_add_reactant_class(
      const ReactionIdBitsets& reactions_bitset_per_reactant, const bool target_only);
  reactant_class_id_t compute_reactant_class_for_species(const species_id_t species_id);

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

  UnimolRxnClassesMap unimol_rxn_class_map;

  BimolRxnClassesMap bimol_rxn_class_map;

  // this map allows to search for reaction classes with ignoring
  // orientation and compartment
  //BimolRxnClassesMap bimol_rxn_class_any_orient_compartment_map;

  reactant_class_id_t next_reactant_class_id;

  // owns reactant classes, indexed by ID
  ReactantClassesVector reactant_classes_vector;
  // contains pointers to objects owned by reactant_classes_vector
  // used to quickly find out whether we already have this reactant class
  std::set<ReactantClass*, ReactantClassLessPtr> reactant_classes_set;

  // indexed by reactant_class_id_t
  std::vector<ReactantClassIdSet> reacting_classes;

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
