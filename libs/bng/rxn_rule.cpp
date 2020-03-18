/*
 * rule.cpp
 *
 *  Created on: Jan 9, 2020
 *      Author: ahusar
 */

#include <iostream>

#include "rxn_rule.h"

using namespace std;

namespace BNG {


void RxnRule::dump_complex_instance_vector(const BNGData& bng_data, const ComplexInstanceVector& complexes) const {

  for (size_t i = 0; i < complexes.size(); i++) {
    complexes[i].dump(bng_data);

    if (i != complexes.size() - 1) {
      cout << " + ";
    }
  }
}


void RxnRule::dump(const BNGData& bng_data) const {
  if (name != "") {
    cout << name << " ";
  }
  dump_complex_instance_vector(bng_data, reactants);

  cout << " -> ";
  dump_complex_instance_vector(bng_data, products);

  cout << " " << rxn_rate;
}

} /* namespace BNG */
