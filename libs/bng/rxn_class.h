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

class SpeciesContainer;

/**
 * Reaction class contains all applicable reactions for a pair of reactants.
 *
 * Reaction classes are created on-the-fly by RxnContainer.
 */
class RxnClass {
public:
  // keeping just IDs, with IDs unlike with pointers we are able to check that the species was 'discarded'
  std::vector<species_id_t> reactants;

  // reactions are owned by RxnContainer
  // order in this vector is important
  std::vector<RxnRule*> reactions;

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
      const SpeciesContainer& all_species_, const BNGConfig& bng_config_,
      const species_id_t reactant_id1, const species_id_t reactant_id2 = SPECIES_ID_INVALID)
    : type(RxnType::Invalid), max_fixed_p(FLT_INVALID), min_noreaction_p(FLT_INVALID),
      all_species(all_species_), bng_config(bng_config_)
    {
    assert(reactant_id1 != SPECIES_ID_INVALID);
    reactants.push_back(reactant_id1);
    if (reactant_id2 != SPECIES_ID_INVALID) {
      reactants.push_back(reactant_id2);
    }
  }

  uint get_num_reactions() const {
    return reactions.size();
  }

  RxnRule* get_rxn(const rxn_index_t rx_index) {
    assert(rx_index < (int)reactions.size());
    return reactions[rx_index];
  }

  orientation_t get_reactant_orientation(uint reactant_index) const {
    assert(!reactions.empty());
    assert(reactant_index < reactions[0]->reactants.size());
    return reactions[0]->reactants[reactant_index].get_orientation();
  }

  void add_rxn_rule(RxnRule* r) {

    // check that the rule was not added already,
    // for now simple pointer comparison
    for (const RxnRule* rxn: reactions) {
      if (r == rxn) {
        // reaction is already present
        return;
      }
    }

    reactions.push_back(r);

    // remember bidirectional mapping for rxn rate updates
    r->add_rxn_class_where_used(this);

    compute_initial_rxn_rates();
  }

  void update_rxn_rates_if_needed(const float_t current_time) {
    // check if any of the reactions needs update
    for (const RxnRule* rxn: reactions) {
      if (rxn->may_update_rxn_rate()) {
        update_variable_rxn_rates(current_time);
        break;
      }
    }
  }

  // this function expects that update_rxn_rates_if_needed was called
  // already for the current time
  float_t get_next_time_of_rxn_rate_update() const {
    float_t min = TIME_FOREVER;
    for (const RxnRule* rxn: reactions) {

      float_t t = rxn->get_next_time_of_rxn_rate_update();
      if (t < min) {
        min = t;
      }
    }
    return min;
  }

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
    return reactants.size() == 1;
  }

  bool is_bimol() const {
    return reactants.size() == 2;
  }

  bool is_absorb_region_border() const {
    return type == RxnType::AbsorbRegionBorder;
  }

  static void dump_array(const std::vector<RxnClass>& vec);

  void dump(const std::string ind = "") const;

private:
  void compute_initial_rxn_rates();

  void update_variable_rxn_rates(const float_t current_time);

  // ----------- MCell-specific -----------
  float_t get_reactant_diffusion(const uint reactant_index) const;
  float_t get_reactant_space_step(const uint reactant_index) const;
  float_t get_reactant_time_step(const uint reactant_index) const;

  float_t compute_pb_factor() const;
  // ^^^^^^^^^^ MCell-specific ^^^^^^^^^^

  // owned by BNGEngine
  const SpeciesContainer& all_species;

  // owned by the simulation engine
  const BNGConfig& bng_config;
};


typedef small_vector<RxnClass*> RxnClassesVector;

} /* namespace BNG */

#endif /* LIBS_BNG_RXN_CLASS_H_ */
