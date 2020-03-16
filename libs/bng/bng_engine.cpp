/*
 * bng_engine.cpp
 *
 *  Created on: Jan 9, 2020
 *      Author: ahusar
 */
#include <iostream>

#include "bng_engine.h"

using namespace std;

namespace BNG {

state_id_t BNGData::find_or_add_state_name(const std::string& s) {
  // rather inefficient search but most probably sufficient for now
  for (state_id_t i = 0; i < state_names.size(); i++) {
    if (state_names[i] == s) {
      return i;
    }
  }

  // not found
  state_names.push_back(s);
  return state_names.size() - 1;
}

component_type_id_t BNGData::find_or_add_component_type(const ComponentType& ct) {
  for (component_type_id_t i = 0; i < component_types.size(); i++) {
    if (component_types[i] == ct) {
      return i;
    }
  }

  // not found
  component_types.push_back(ct);
  return component_types.size() - 1;
}

molecule_type_id_t BNGData::find_or_add_molecule_type(const MoleculeType& mt) {
  for (component_type_id_t i = 0; i < molecule_types.size(); i++) {
    if (molecule_types[i] == mt) {
      return i;
    }
  }

  // not found
  molecule_types.push_back(mt);
  return molecule_types.size() - 1;
}


void BNGData::dump(const bool as_bngl) {
  if (as_bngl) {
    cout << "begin molecule types\n";

    for (const MoleculeType& mt: molecule_types) {

      cout << "  " << mt.name << "(";

      for (size_t i = 0; i < mt.component_type_ids.size(); i++) {

        const ComponentType& ct = get_component_type(mt.component_type_ids[i]);

        cout << ct.name;

        for (state_id_t state_id: ct.allowed_state_ids) {
          cout << "~" << get_state_name(state_id);
        }
        if (i != mt.component_type_ids.size() - 1) {
          cout << ", ";
        }
      }

      cout << ")\n";
    }

    cout << "end molecule types\n";
  }
  else {
    assert(false && "TODO");
  }
}


} /* namespace BNG */
