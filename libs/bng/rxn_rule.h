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

  // the complex species are patterns
  //
  // there is a potential for optimizations, e.g.
  // to make a set of species that match the patterns, but let's keep it
  // for later
  ComplexSpeciesInstanceVector reactants;
  ComplexSpeciesInstanceVector products;

  float_t reaction_rate;

  bool operator ==(const RxnRule& rr2) {
    // ordering of components in a molecule is important
    // two component types must have the same id, this is ensured in find_or_add_component_type
    return
        name == rr2.name &&
        reactants == rr2.reactants && products == rr2.products &&
        reaction_rate == rr2.reaction_rate;
  }
};

} /* namespace BNG */

#endif /* LIBS_BNG_RXN_RULE_H_ */
