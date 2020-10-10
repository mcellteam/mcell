/*
 * bng_data.cpp
 *
 *  Created on: Mar 24, 2020
 *      Author: ahusar
 */

#include <iostream>

#include "bng/bng_data.h"

using namespace std;

namespace BNG {

void BNGData::clear() {
  state_names.clear();
  component_types.clear();
  molecule_types.clear();
  rxn_rules.clear();
}


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


// may return STATE_ID_INVALID when the name was not found
state_id_t BNGData::find_state_id(const std::string& name) const {
  for (state_id_t i = 0; i < state_names.size(); i++) {
    if (state_names[i] == name) {
      return i;
    }
  }
  return MOL_TYPE_ID_INVALID;
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


component_type_id_t BNGData::find_component_type_id(const MolType& mt, const std::string& name) const {
  for (component_type_id_t ct_id: mt.component_type_ids) {
    if (component_types[ct_id].name == name) {
      return ct_id;
    }
  }
  return COMPONENT_TYPE_ID_INVALID;
}


// compares only name
mol_type_id_t BNGData::find_or_add_molecule_type(const MolType& mt) {
  // TODO LATER: check that if there is a molecule type with the same name,
  // it a subset of the same components and allowed states

  mol_type_id_t existing_mt_id = find_molecule_type_id(mt.name);
  if (existing_mt_id != MOL_TYPE_ID_INVALID) {
    return existing_mt_id;
  }

  // not found
  molecule_types.push_back(mt);
  molecule_types.back().set_finalized();
  return molecule_types.size() - 1;
}


// may return MOLECULE_TYPE_ID_INVALID when the name was not found
mol_type_id_t BNGData::find_molecule_type_id(const std::string& name) const {
  for (mol_type_id_t i = 0; i < molecule_types.size(); i++) {
    const MolType& mt = molecule_types[i];
    if (mt.name == name) {
      return i;
    }
  }
  return MOL_TYPE_ID_INVALID;
}


// asserts if compartment with the same name already exists
compartment_id_t BNGData::add_compartment(const Compartment& c) {
  assert(find_compartment_id(c.name) == COMPARTMENT_ID_INVALID);
  compartment_id_t id = compartments.size();
  compartments.push_back(c);
  compartments.back().id = id;
  return id;
}


// returns COMPARTMENT_ID_INVALID if compartment with this name was not found
compartment_id_t BNGData::find_compartment_id(const std::string& name) const {
  for (compartment_id_t i = 0; i < compartments.size(); i++) {
    if (compartments[i].name == name) {
      return i;
    }
  }
  return COMPARTMENT_ID_INVALID;
}


Compartment* BNGData::find_compartment(const std::string& name) {
  compartment_id_t id = find_compartment_id(name);
  if (id == COMPARTMENT_ID_INVALID) {
    return nullptr;
  }
  else {
    return &get_compartment(id);
  }
}

const Compartment* BNGData::find_compartment(const std::string& name) const {
  compartment_id_t id = find_compartment_id(name);
  if (id == COMPARTMENT_ID_INVALID) {
    return nullptr;
  }
  else {
    return &get_compartment(id);
  }
}


rxn_rule_id_t BNGData::find_or_add_rxn_rule(const RxnRule& rr) {
  // TODO LATER: check that if there is a reaction with the same
  //       reactants and products that the reaction rate is the same
  for (rxn_rule_id_t i = 0; i < rxn_rules.size(); i++) {
    if (rxn_rules[i] == rr) {
      return i;
    }
  }

  // not found
  rxn_rule_id_t id = rxn_rules.size();
  rxn_rules.push_back(rr);
  rxn_rules.back().id = id;
  return id;
}


void BNGData::dump_molecule_types_as_bngl() {
  cout << "begin molecule types\n";

  for (const MolType& mt: molecule_types) {
    cout << "  ";
    mt.dump(*this);
    cout << "\n";
  }

  cout << "end molecule types\n";
}


void BNGData::dump_reaction_rules_as_bngl() {
  cout << "begin reaction rules\n";

  for (const RxnRule& rr: rxn_rules) {
    cout << "  ";
    rr.dump();
    cout << "\n";
  }

  cout << "end reaction rules\n";

}


void BNGData::dump() {
  dump_molecule_types_as_bngl();
  dump_reaction_rules_as_bngl();
}

} /* namespace BNG2 */
