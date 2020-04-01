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

#include "bng_defines.h"

#include "rxn_rule.h"

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
  RxnClass(const SpeciesContainer& all_species_, const species_id_t reactant_id1, const species_id_t reactant_id2 = SPECIES_ID_INVALID)
    : type(RxnType::Invalid), max_fixed_p(FLT_INVALID), min_noreaction_p(FLT_INVALID), all_species(all_species_)
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

  const RxnRule* get_reaction(const reaction_index_t rx_index) const {
    assert(rx_index < (int)reactions.size());
    return reactions[rx_index];
  }

  orientation_t get_reactant_orientation(uint reactant_index) const {
    assert(!reactions.empty());
    assert(reactant_index < reactions[0]->reactants.size());
    return reactions[0]->reactants[reactant_index].get_orientation();
  }

  void add_rxn_rule(const BNGConfig& bng_config, RxnRule* r) {

    // check that the rule was not added already,
    // for now simple pointer comparison
    for (const RxnRule* rxn: reactions) {
      if (r == rxn) {
        // reaction is already present
        return;
      }
    }

    reactions.push_back(r);

    update(bng_config);
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

  bool is_absorb() const {
    return type == RxnType::AbsorbRegionBorder;
  }

  static void dump_array(const BNGData& bng_data, const std::vector<RxnClass>& vec);

  void dump(const BNGData& bng_data, const std::string ind = "") const;

private:
  void update(const BNGConfig& bng_config);

  // ----------- MCell-specific -----------
  float_t get_reactant_diffusion(const uint reactant_index) const;
  float_t get_reactant_space_step(const uint reactant_index) const;
  float_t get_reactant_time_step(const uint reactant_index) const;

  float_t compute_pb_factor(const BNGConfig& bng_config) const;
  // ^^^^^^^^^^ MCell-specific ^^^^^^^^^^

  // owned by BNGEngine
  const SpeciesContainer& all_species;
};


typedef small_vector<const RxnClass*> RxnClassesVector;

} /* namespace BNG */

#endif /* LIBS_BNG_RXN_CLASS_H_ */
