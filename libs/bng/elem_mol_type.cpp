/*
 * elementary_molecule.cpp
 *
 *  Created on: Jan 9, 2020
 *      Author: ahusar
 */
#include <bng/elem_mol_type.h>
#include <iostream>
#include <sstream>

#include "bng/bng_engine.h"

using namespace std;

namespace BNG {

// ------------- ComponentType -------------
std::string ComponentType::to_str(const BNGData& bng_data) const {
  stringstream out;
  out << name;

  for (state_id_t state_id: allowed_state_ids) {
    out << "~" << bng_data.get_state_name(state_id);
  }
  return out.str();
}


void ComponentType::dump(const BNGData& bng_data) const {
  cout << to_str(bng_data);
}


// ------------- MoleculeType -------------
std::string ElemMolType::to_str(const BNGData& bng_data) const {
  stringstream out;

  out << name << "(";

  for (size_t i = 0; i < component_type_ids.size(); i++) {

    const ComponentType& ct = bng_data.get_component_type(component_type_ids[i]);
    out << ct.to_str(bng_data);

    if (i != component_type_ids.size() - 1) {
      out << ",";
    }
  }

  out << ")";
  return out.str();
}


void ElemMolType::dump(const BNGData& bng_data) const {

  cout << to_str(bng_data);
  cout <<
      " D=" << D <<
      ", custom_time_step=" << custom_time_step <<
      ", custom_space_step=" << custom_space_step;
}

} /* namespace BNG */
