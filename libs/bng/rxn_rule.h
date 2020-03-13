/*
 * rule.h
 *
 *  Created on: Jan 9, 2020
 *      Author: ahusar
 */

#ifndef LIBS_BNG_RXN_RULE_H_
#define LIBS_BNG_RXN_RULE_H_

#include <string>

#include "bng_defines.h"

#include "complex_species.h"

namespace BNG {

// BNG reaction rule
// rules are only unidirectional,
// if there is a reversible reaction in BNGL definition,
// two
class RxnRule {
public:
  std::string name;

  small_vector<ComplexSpecies> reactants;
  small_vector<ComplexSpecies> products;

  float_t reaction_rate;

};

} /* namespace BNG */

#endif /* LIBS_BNG_RXN_RULE_H_ */
