/*
 * RxnClass.cpp
 *
 *  Created on: Mar 25, 2020
 *      Author: ahusar
 */

#include <iostream>

#include "rxn_class.h"

using namespace std;

namespace BNG {


void RxnClass::dump_array(const BNGData& bng_data, const vector<RxnClass>& vec) {
  cout << "Reaction class array: " << (vec.empty() ? "EMPTY" : "") << "\n";

  for (size_t i = 0; i < vec.size(); i++) {
    cout << i << ":\n";
    vec[i].dump(bng_data, "  ");
  }
}

void RxnClass::dump(const BNGData& bng_data, const std::string ind) const {
  cout << ind << "max_fixed_p: \t\t" << max_fixed_p << " [float_t] \t\t\n";
  cout << ind << "min_noreaction_p: \t\t" << min_noreaction_p << " [float_t] \t\t\n";

  for (const RxnRule* rxn: reactions) {
    rxn->dump(bng_data, ind);
  }
}

} /* namespace BNG */
