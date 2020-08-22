/*
 * RxnClass.h
 *
 *  Created on: Mar 25, 2020
 *      Author: ahusar
 */

#ifndef LIBS_BNG_RXN_CLASS_H_
#define LIBS_BNG_RXN_CLASS_H_

#include <string>
#include <iostream>

#include "bng/bng_defines.h"
#include "bng/rxn_rule.h"

namespace BNG {

class RxnContainer;
class SpeciesContainer;

/**
 * Used to hold information on the probability and products for a given reaction rule.
 * One rxn rule can possibly have more different products,
 * e.g. when applying rule A(b~0) -> A(b~1) on reactant A(a~0,b~0).A(a~1,b~0),
 * the result can be either A(a~0,b~1).A(a~1,b~0) or A(a~0,b~0).A(a~1,b~1).
 * So for each of the possible products, one RxnPathway is created in a RxnClass.
 */
class RxnClassPathway {
public:
  // ID of rxn rule from which this pathway was created
  rxn_rule_id_t rxn_rule_id;

  // cumulative probability
  float_t cum_prob;

  // probability for this specific pathway to be selected when
  // reactants interact (or when an unimol rxn is executed)
  float_t pathway_prob;

  // specific variant of products for this pathway
  std::vector<species_id_t> product_species;
};

/**
 * Reaction class contains all applicable reactions for a pair of reactants.
 *
 * Reaction classes are created on-the-fly by RxnContainer.
 */
class RxnClass {
public:
  // keeping just IDs, with IDs unlike with pointers we are able to check that the species was 'discarded'
  // the rxn class was created specifically for a single or a pair of reactants, therefore
  // the species id is is a specific instance of a complex that matches our pattern
  // and if there are multiple species that match, multiple rxn classes are created
  std::vector<species_id_t> specific_reactants;

private:
  // reactions are owned by RxnContainer, order in this vector is important
  std::vector<rxn_rule_id_t> rxn_rule_ids;

  // initialized in ...?
  std::vector<RxnClassPathway> pathways;

public:
  // Standard reaction or special such as Reflect, Transparent or Absorb
  RxnType type;

  // ----------- MCell-specific -----------
  // Maximum 'p' for region of p-space for all non-cooperative pathways
  float_t max_fixed_p;

  // Minimum 'p' for region of p-space which is always in the non-reacting "pathway". (note that
  // cooperativity may mean that some values of p less than this still do not produce a reaction)
  float_t min_noreaction_p;

  // Cumulative probabilities for specific reactions, based on all reactions of the class
  // has same size as reactions
  std::vector<float_t> cum_probs;
  // ^^^^^^^^^^ MCell-specific ^^^^^^^^^^

public:
  RxnClass(
      RxnContainer& all_rxns_, const SpeciesContainer& all_species_, const BNGConfig& bng_config_,
      const species_id_t reactant_id1, const species_id_t reactant_id2 = SPECIES_ID_INVALID)
    : type(RxnType::Invalid), max_fixed_p(FLT_INVALID), min_noreaction_p(FLT_INVALID),
      all_rxns(all_rxns_), all_species(all_species_), bng_config(bng_config_), bimol_vol_rxn_flag(false)
    {
    assert(reactant_id1 != SPECIES_ID_INVALID);
    specific_reactants.push_back(reactant_id1);
    if (reactant_id2 != SPECIES_ID_INVALID) {
      specific_reactants.push_back(reactant_id2);
    }
  }

  uint get_num_reactions() const {
    return rxn_rule_ids.size();
  }

  RxnRule* get_rxn_for_pathway(const rxn_class_pathway_index_t pathway_index);

  /*const std::vector<RxnRule*>& get_rxns() const {
    return rxn_rules;
  }*/

  orientation_t get_reactant_orientation(uint reactant_index) const;

  // does not do pathways update
  void add_rxn_rule_no_update(RxnRule* r) {

    if (rxn_rule_ids.empty()) {
      bimol_vol_rxn_flag = r->is_bimol_vol_rxn();
    }
    else {
      assert(bimol_vol_rxn_flag == r->is_bimol_vol_rxn());
    }

    // check that the rule was not added already,
    // for now simple pointer comparison
    for (rxn_rule_id_t id: rxn_rule_ids) {
      if (r->id == id) {
        // reaction is already present
        return;
      }
    }

    rxn_rule_ids.push_back(r->id);

    // remember bidirectional mapping for rxn rate updates
    r->add_rxn_class_where_used(this);
  }

  // must be called after all rxn rules were added to this reaction
  // class called automatically from update_rxn_rates_if_needed
  void update_rxn_pathways();

  void update_rxn_rates_if_needed(const float_t current_time);

  // this function expects that update_rxn_rates_if_needed was called
  // already for the current time
  float_t get_next_time_of_rxn_rate_update() const;

  bool is_standard() const {
    return type == RxnType::Standard;
  }

  // there is exactly one reaction for this type
  bool is_reflect() const {
    return type == RxnType::Reflect;
  }

  bool is_transparent() const {
    return type == RxnType::Transparent;
  }

  bool is_unimol() const {
    return specific_reactants.size() == 1;
  }

  bool is_bimol() const {
    return specific_reactants.size() == 2;
  }

  bool is_bimol_vol_rxn_class() const {
    #ifndef NDEBUG
      debug_check_bimol_vol_rxn_flag();
    #endif
    return bimol_vol_rxn_flag;
  }

  bool is_absorb_region_border() const {
    return type == RxnType::AbsorbRegionBorder;
  }

  // mostly for debug purposes
  bool is_simple() const;

  species_id_t get_second_species_id(species_id_t reactant_id) const {
    assert(is_bimol());
    if (specific_reactants[0] != reactant_id) {
      return specific_reactants[0];
    }
    else {
      return specific_reactants[1];
    }
  }

  static void dump_array(const std::vector<RxnClass>& vec);

  void dump(const std::string ind = "") const;

private:
  void update_variable_rxn_rates(const float_t current_time);

  void debug_check_bimol_vol_rxn_flag() const;

  // ----------- MCell-specific -----------
  float_t get_reactant_diffusion(const uint reactant_index) const;
  float_t get_reactant_space_step(const uint reactant_index) const;
  float_t get_reactant_time_step(const uint reactant_index) const;

  float_t compute_pb_factor() const;
  // ^^^^^^^^^^ MCell-specific ^^^^^^^^^^

  // owned by BNGEngine
  RxnContainer& all_rxns;
  const SpeciesContainer& all_species;

  // owned by the simulation engine
  const BNGConfig& bng_config;

  // flag for optimized testing of this rxn class
  bool bimol_vol_rxn_flag;
};


typedef small_vector<RxnClass*> RxnClassesVector;

} /* namespace BNG */

#endif /* LIBS_BNG_RXN_CLASS_H_ */
