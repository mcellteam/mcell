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

class BNGData;

struct ComplexMoleculeIndex {
  uint complex_index;
  uint molecule_index;
};

struct ComplexMoleculeInstancePair {
  ComplexMoleculeIndex cmi1;
  ComplexMoleculeIndex cmi2;
};

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
  ComplexInstanceVector reactants;
  ComplexInstanceVector products;

  // set to true if it was possible to do a mapping between reactants and products
  bool molecule_instances_are_maintained;

  // matching between molecules of reactants and molecules of products,
  // contains information only if molecule_instances_are_maintained is true
  small_vector<ComplexMoleculeInstancePair> mapping;


  float_t rxn_rate;

  bool operator ==(const RxnRule& rr2) {
    // ordering of components in a molecule is important
    // two component types must have the same id, this is ensured in find_or_add_component_type
    return
        name == rr2.name &&
        reactants == rr2.reactants && products == rr2.products &&
        rxn_rate == rr2.rxn_rate;
  }

  void dump(const BNGData& bng_data) const;

  // checks if it is possible to create a mapping from reactants to products and
  // sets members molecule_instances_are_maintained and mapping,
  // might write some error messages to the msgs stream,
  // returns true if errors were encountered
  bool compute_reactants_products_mapping(std::stringstream& msgs);

private:
  void dump_complex_instance_vector(
      const BNGData& bng_data,
      const ComplexInstanceVector& complexes) const;
};

} /* namespace BNG */

#endif /* LIBS_BNG_RXN_RULE_H_ */
