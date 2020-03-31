/*
 * bng_data.h
 *
 *  Created on: Mar 24, 2020
 *      Author: ahusar
 */

#ifndef LIBS_BNG_BNG_DATA_H_
#define LIBS_BNG_BNG_DATA_H_

#include "bng_defines.h"
#include "cplx_instance.h"

#include "mol_type.h"
#include "rxn_rule.h"


namespace BNG {

/**
 * Data shared among all instances of BNGEngines
 * Usually constant, initialized when BNGL is parsed
 */
class BNGData {
private:
  // indexed with state_id_t
  std::vector<std::string> state_names;

  // indexed with component_type_id_t
  std::vector<ComponentType> component_types;

  // indexed with molecule_type_id_t
  std::vector<MolType> molecule_types;

  // indexed with rxn_rule_id_t
  std::vector<RxnRule> rxn_rules;
  
public:
  // -------- component state names --------

  state_id_t find_or_add_state_name(const std::string& s);

  // may return STATE_ID_INVALID when the name was not found
  state_id_t find_state_id(const std::string& name) const;

  const std::string& get_state_name(const state_id_t id) const {
    assert(id < state_names.size());
    return state_names[id];
  }


  // -------- components --------

  component_type_id_t find_or_add_component_type(const ComponentType& ct);

  const ComponentType& get_component_type(const component_type_id_t id) const {
    assert(id < component_types.size());
    return component_types[id];
  }


  // -------- molecule types --------

  mol_type_id_t find_or_add_molecule_type(const MolType& mt);

  // may return MOLECULE_TYPE_ID_INVALID when the name was not found
  mol_type_id_t find_molecule_type_id(const std::string& name) const;

  const MolType& get_molecule_type(const mol_type_id_t id) const {
    assert(id < molecule_types.size());
    return molecule_types[id];
  }


  // -------- reaction rules --------

  rxn_rule_id_t find_or_add_rxn_rule(const RxnRule& rr);

  // -------- utilities --------
  void dump();

private:
  void dump_molecule_types_as_bngl();
  void dump_reaction_rules_as_bngl();
};



} /* namespace BNG2 */

#endif /* LIBS_BNG_BNG_DATA_H_ */
