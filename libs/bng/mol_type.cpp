/*
 * elementary_molecule.cpp
 *
 *  Created on: Jan 9, 2020
 *      Author: ahusar
 */
#include "bng/mol_type.h"

#include <iostream>

#include "bng/bng_engine.h"

using namespace std;

namespace BNG {

// ------------- ComponentType -------------
void ComponentType::dump(const BNGData& bng_data) const {
  cout << name;

  for (state_id_t state_id: allowed_state_ids) {
    cout << "~" << bng_data.get_state_name(state_id);
  }
}


// ------------- MoleculeType -------------
void MolType::dump(const BNGData& bng_data) const {
  cout << name << "(";

  for (size_t i = 0; i < component_type_ids.size(); i++) {

    const ComponentType& ct = bng_data.get_component_type(component_type_ids[i]);
    ct.dump(bng_data);

    if (i != component_type_ids.size() - 1) {
      cout << ", ";
    }
  }

  cout << ")";
  cout <<
      ", D=" << D <<
      ", custom_time_step=" << custom_time_step <<
      ", custom_space_step=" << custom_space_step;
}

} /* namespace BNG */
