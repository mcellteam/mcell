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

enum class RxnClassType {
  Standard, // any other reaction than below
  Transparent,
  Reflect,
  AbsorbRegionBorder
};


// corresponds to reaction in mcell3 code
// coupled reactions for given species, serves for caching
// created on-the-fly
class RxnClass {
public:
  RxnClassType type;

  /* Maximum 'p' for region of p-space for all non-cooperative pathways */
  float_t max_fixed_p;

  /* Minimum 'p' for region of p-space which is always in the non-reacting "pathway". (note that
     cooperativity may mean that some values of p less than this still do not produce a reaction) */
  float_t min_noreaction_p;

  // reactants are copied into specific reactions as well
  // because a different order might be needed
  CplxInstanceVector reactants;

  // reactions are owned by RxnContainer
  std::vector<RxnRule*> reactions;

  // Cumulative probabilities for specific reactions, based on all reactions of the class
  // has same size as reactions
  std::vector<float_t> cum_probs;

public:

  uint get_num_reactions() const {
    return reactions.size();
  }

  const RxnRule* get_reaction(reaction_index_t rx_index) const {
    assert(rx_index < (int)reactions.size());
    return reactions[rx_index];
  }

  void add_and_finalize_reaction(RxnRule* r, const float_t cum_prob) {
    assert(reactions.size() == cum_probs.size());

    reactions.push_back(r);
    reactions.back()->finalize();

    cum_probs.push_back(cum_prob);
  }

  // there are no pathways for this type of reactions
  bool is_reflect() const {
    return type == RxnClassType::Reflect;
  }

  bool is_transparent() const {
    return type == RxnClassType::Transparent;
  }

  bool is_absorb() const {
    return type == RxnClassType::AbsorbRegionBorder;
  }

  static void dump_array(const BNGData& bng_data, const std::vector<RxnClass>& vec);

  void dump(const BNGData& bng_data, const std::string ind = "") const;
};


typedef small_vector<const RxnClass*> RxnClassesVector;

} /* namespace BNG */

#endif /* LIBS_BNG_RXN_CLASS_H_ */
