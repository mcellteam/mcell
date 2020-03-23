/*
 * bng_engine.h
 *
 *  Created on: Jan 9, 2020
 *      Author: ahusar
 */

#ifndef LIBS_BNG_BNG_ENGINE_H_
#define LIBS_BNG_BNG_ENGINE_H_

#include "bng_defines.h"
#include "cplx_instance.h"

#include "mol_type.h"
#include "rxn_rule.h"

namespace BNG {

// TODO: split, make bng_data.h file
// Data shared among all instances of BNGEngines
// Usually constant, initialized when BNGL is parsed
// TODO: make arrays private
class BNGData {
public:
  // -------- component states --------

  // indexed with state_id_t
  std::vector<std::string> state_names;

  state_id_t find_or_add_state_name(const std::string& s);

  // may return STATE_ID_INVALID when the name was not found
  state_id_t find_state_id(const std::string& name) const;

  const std::string& get_state_name(const state_id_t id) const {
    assert(id < state_names.size());
    return state_names[id];
  }


  // -------- components --------

  // indexed with component_type_id_t
  std::vector<ComponentType> component_types;

  component_type_id_t find_or_add_component_type(const ComponentType& ct);

  const ComponentType& get_component_type(const component_type_id_t id) const {
    assert(id < component_types.size());
    return component_types[id];
  }


  // -------- molecule types --------

  // indexed with molecule_type_id_t
  std::vector<MolType> molecule_types;

  mol_type_id_t find_or_add_molecule_type(const MolType& mt);

  // may return MOLECULE_TYPE_ID_INVALID when the name was not found
  mol_type_id_t find_molecule_type_id(const std::string& name) const;

  const MolType& get_molecule_type(const mol_type_id_t id) const {
    assert(id < molecule_types.size());
    return molecule_types[id];
  }

  

  // -------- reaction rules --------
  std::vector<RxnRule> rxn_rules;

  rxn_rule_id_t find_or_add_rxn_rule(const RxnRule& rr);


  // -------- utilities --------
  void dump(const bool as_bngl = true);

private:
  void dump_molecule_types_as_bngl();
  void dump_reaction_rules_as_bngl();
};



// Each partition will have its own BNG engine,
// the contents might change a lot during execution and string comparison
// on using BNGL components we still have a way how to unify all the instances
class BNGEngine {
public:


  //component_index_t next_component_index;

  // checks cache - 2 types of caches can/cannot react
  // if not found, we need to decide
  void can_react(
      cplx_species_id_t a,
      cplx_species_id_t b,
      small_vector<cplx_species_id_t>& reactions
  );

  void get_reactants(
      cplx_species_id_t b,
      small_vector<cplx_species_id_t>& reactions
  );

  // search whether two molecules can react is done
  //std::vector<ComplexSpecies> complex_species;


  // cache of complex species indices that can interact together

  const BNGData& bng_data;
};



} /* namespace BNG */

#endif /* LIBS_BNG_BNG_ENGINE_H_ */
