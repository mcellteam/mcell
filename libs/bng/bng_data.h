/*
 * bng_data.h
 *
 *  Created on: Mar 24, 2020
 *      Author: ahusar
 */

#ifndef LIBS_BNG_BNG_DATA_H_
#define LIBS_BNG_BNG_DATA_H_

#include "bng/bng_defines.h"
#include "bng/mol_type.h"
#include "bng/rxn_rule.h"
#include "bng/cplx.h"

namespace BNG {


class SeedSpecies {
public:
  SeedSpecies(const BNGData* bng_data)
    : cplx(bng_data),
      count(0),
      compartment_id(COMPARTMENT_ID_INVALID) {
  }

  Cplx cplx;
  float_t count; // count of molecules to be released (exact value from the BNGL file)
  compartment_id_t compartment_id; // compartment where molecules should be released
};


enum class ObservableType {
  Invalid,
  Species,
  Molecules
};


class Observable {
public:
  Observable()
    : type(ObservableType::Invalid) {
  }

  ObservableType type;
  std::string name;
  std::vector<Cplx> patterns;
};


class Compartment {
public:
  Compartment()
    : is_3d(true), volume(FLT_INVALID) {
  }
  std::string name;
  bool is_3d; // 2d if this member is false
  float_t volume;
  std::string parent_name;
};


/**
 * Data shared among all instances of BNGEngines
 * Usually constant, initialized when BNGL is parsed
 */
class BNGData {
private:
  // indexed with state_id_t
  std::vector<std::string> state_names;

  // indexed with component_type_id_t,
  // contains unique component types even though their name may be identical
  std::vector<ComponentType> component_types;

  // indexed with molecule_type_id_t
  std::vector<MolType> molecule_types;

  // indexed with compartment_id_t
  std::vector<Compartment> compartments;

  // indexed with rxn_rule_id_t
  // rxn rules are then in rxn container, this is a temporary placeholder
  // for parsing result
  std::vector<RxnRule> rxn_rules;
  
  // not referenced by any other data in this class,
  // keeping for cases when the parameter values might be needed for other purposes
  std::map<std::string, float_t> parameters;

  // contents of the seed species section
  // not used directly but can be converted to other representations
  std::vector<SeedSpecies> seed_species;

  // contents of the observables section
  // not used directly but can be converted to other representations
  std::vector<Observable> observables;

public:
  void clear();

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

  // may return COMPONENT_TYPE_ID_INVALID when the name was not found
  // among components allowed for this molecule type
  component_type_id_t find_component_type_id(const MolType& mt, const std::string& name) const;

  const ComponentType& get_component_type(const component_type_id_t id) const {
    assert(id < component_types.size());
    return component_types[id];
  }

  ComponentType& get_component_type(const component_type_id_t id) {
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

  MolType& get_molecule_type(const mol_type_id_t id) {
    assert(id < molecule_types.size());
    return molecule_types[id];
  }

  const std::vector<MolType>& get_molecule_types() const {
    return molecule_types;
  }

  // -------- compartments --------
  // asserts if compartment with the same name already exists
  compartment_id_t add_compartment(const Compartment& c);

  // returns COMPARTMENT_ID_INVALID if compartment with this name was not found
  compartment_id_t find_compartment_id(const std::string& name) const;

  // returns nullptr if compartment with this name was not found
  const Compartment* find_compartment(const std::string& name) const;

  const Compartment& get_compartment(const compartment_id_t id) const {
    assert(id < compartments.size());
    return compartments[id];
  }

  const std::vector<Compartment>& get_compartments() const {
    return compartments;
  }
  // -------- reaction rules --------

  rxn_rule_id_t find_or_add_rxn_rule(const RxnRule& rr);

  const std::vector<RxnRule>& get_rxn_rules() const {
    return rxn_rules;
  }


  // -------- seed species --------

  void add_seed_species(const SeedSpecies& ss) {
    seed_species.push_back(ss);
  }

  const std::vector<SeedSpecies>& get_seed_species() const {
    return seed_species;
  }

  // -------- observables --------

  void add_observable(const Observable& o) {
    observables.push_back(o);
  }

  const std::vector<Observable>& get_observables() const {
    return observables;
  }

  // -------- parameters --------

  void add_parameter(const std::string& name, const float_t& value) {
    assert(parameters.count(name) == 0);
    parameters[name] = value;
  }

  // returns true and sets value if found, returns falser otherwise
  bool get_parameter_value(const std::string& name, float_t& value) const {
    auto it = parameters.find(name);
    if (it == parameters.end()) {
      return false;
    }
    else {
      value = it->second;
      return true;
    }
  }

  const std::map<std::string, float_t>& get_parameters() const {
    return parameters;
  }

  // -------- utilities --------
  void dump();

private:
  void dump_molecule_types_as_bngl();
  void dump_reaction_rules_as_bngl();
};

} /* namespace BNG2 */

#endif /* LIBS_BNG_BNG_DATA_H_ */
