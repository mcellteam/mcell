/*
 * complex_species.cpp
 *
 *  Created on: Jan 9, 2020
 *      Author: ahusar
 */
#include <iostream>

#include "ast.h"
#include "bng_engine.h"
#include "cplx_instance.h"
#include "mol_type.h"

using namespace std;

namespace BNG {

// ------------- ComponentInstance -------------
void ComponentInstance::dump(const BNGData& bng_data) const {

  const ComponentType& ct = bng_data.get_component_type(component_type_id);
  cout << ct.name;

  assert(state_id != STATE_ID_INVALID);
  if (state_id != STATE_ID_DONT_CARE) {
    cout << "~" << bng_data.get_state_name(state_id);
  }

  assert(state_id != BOND_VALUE_INVALID);
  if (bond_value == BOND_VALUE_ANY) {
    cout << "!" + BOND_STR_ANY;
  }
  else if (bond_value != BOND_VALUE_NO_BOND) {
    cout << "!"  << bond_value;
  }
}


// ------------- MoleculeInstance -------------
void MolInstance::initialize_components_types(const MolType& mt) {
  for (component_type_id_t component_type_id: mt.component_type_ids) {
    // state is don't care, no bond
    component_instances.push_back(component_type_id);
  }
}


// searches for component with name
uint MolInstance::get_corresponding_component_index(
    const BNGData& bng_data,
    const MolType& mt,
    const std::string& name,
    const uint starting_index
) const {
  for (uint i = starting_index; i < mt.component_type_ids.size(); i++) {
    const ComponentType& ct = bng_data.get_component_type(mt.component_type_ids[i]);

    // we are looking for the first component with the same name
    if (ct.name == name) {
      return i;
    }
  }

  return INDEX_INVALID;
}


void MolInstance::dump(const BNGData& bng_data) const {
  const MolType& mt = bng_data.get_molecule_type(mol_type_id);
  cout << mt.name << "(";

  for (size_t i = 0; i < component_instances.size(); i++) {
    component_instances[i].dump(bng_data);
    if (i != component_instances.size() - 1) {
      cout << ", ";
    }
  }
  cout << ")";
}


// ------------- ComplexInstance -------------
void CplxInstance::dump(const BNGData& bng_data) const {
  for (size_t i = 0; i < mol_patterns.size(); i++) {
    mol_patterns[i].dump(bng_data);

    if (i != mol_patterns.size() - 1) {
      cout << ".";
    }
  }
}

} /* namespace BNG */
