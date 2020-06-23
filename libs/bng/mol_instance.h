/*
 * complex_species.h
 *
 *  Created on: Jan 9, 2020
 *      Author: ahusar
 */

#ifndef LIBS_BNG_MOL_INSTANCE_H_
#define LIBS_BNG_MOL_INSTANCE_H_

#include <iostream>

#include "bng/bng_defines.h"
#include "bng/base_flag.h"

namespace BNG {

class MolType;
class BNGData;

/**
 * Pattern vs instance:
 * - patterns are used in reaction rules, they are written in the don't care don't write manner and
 *   may use '+' as bond
 *
 * - instances are fully specified, each component has a given state (if it has states)
 */
class ComponentInstance {
public:
  ComponentInstance()
    : component_type_id(COMPONENT_TYPE_ID_INVALID), state_id(STATE_ID_DONT_CARE), bond_value(BOND_VALUE_ANY),
      explicitly_listed_in_pattern(false) {
  }

  ComponentInstance(const component_type_id_t id)
    : component_type_id(id), state_id(STATE_ID_DONT_CARE), bond_value(BOND_VALUE_ANY),
      explicitly_listed_in_pattern(false) {
  }

  ComponentInstance(const ComponentInstance& other)
    : component_type_id(other.component_type_id),
      state_id(other.state_id),
      bond_value(other.bond_value),
      explicitly_listed_in_pattern(other.explicitly_listed_in_pattern)
    {
  }

  // type id can be also found from parent's MoleculeInstance::molecule_type_id
  // maybe remove
  component_type_id_t component_type_id;

  // specific state or STATE_ID_DONT_CARE if we don't care or unset
  state_id_t state_id;

  // number for bonds (0, 1, ...),
  // BOND_VALUE_ANY for patterns (+),
  // BOND_VALUE_NO_BOND if this component has no bond
  bond_value_t bond_value;

  // true if this component was explicitly listed in a pattern
  bool explicitly_listed_in_pattern;

  bool operator == (const ComponentInstance& comp2) const  {
    return
        component_type_id == comp2.component_type_id &&
        state_id == comp2.state_id &&
        bond_value == comp2.bond_value;
  }

  bool bond_has_numeric_value() const {
    assert(bond_value != BOND_VALUE_INVALID);
    return bond_value != BOND_VALUE_ANY && bond_value != BOND_VALUE_NO_BOND;
  }

  bool state_is_set() const {
    assert(state_id != STATE_ID_INVALID);
    return state_id != STATE_ID_DONT_CARE;
  }

  std::string to_str(const BNGData& bng_data) const;
  void dump(const BNGData& bng_data, const std::string& ind = "") const;
};


/**
 * Instance of a molecule type.
 * Similarly as ComplexSpeciesInstance, it can be used as instance and as pattern.
 */
class MolInstance: public BaseSpeciesCplxMolFlag {
public:
  // ID of this molecule type in BNGData::molecule_types
  mol_type_id_t mol_type_id;

  // has the same number of elements as MoleculeType::component_type_ids
  small_vector<ComponentInstance> component_instances;

public:
  MolInstance()
    : mol_type_id(MOL_TYPE_ID_INVALID) {
  }

  MolInstance(const MolInstance& other)
    : mol_type_id(other.mol_type_id),
      component_instances(other.component_instances) {
  }

  void finalize_flags() {
    set_finalized();
    // flag about molecule type must be set
    assert(is_vol() || is_surf() || is_reactive_surface());
  }

  // returns true if this object as a pattern matches second instance
  bool matches_simple(const MolInstance& inst) const {
    assert(component_instances.size() == 0 && inst.component_instances.size() == 0 &&
        "Method can be used only for simple complexes, i.e. without components.");

    return mol_type_id == inst.mol_type_id;
  }

  bool operator == (const MolInstance& other) const  {
    return
        mol_type_id == other.mol_type_id &&
        component_instances == other.component_instances &&
        get_flags() == other.get_flags();
  }

  void initialize_components_types(const MolType& mt);

  // searches for component with name
  uint get_corresponding_component_index(
      const BNGData& bng_data,
      const MolType& mt,
      const std::string& name,
      const uint starting_index
  ) const;

  std::string to_str(const BNGData& bng_data, const bool only_explicit = false) const;
  //void dump(const BNGData& bng_data, const bool for_diff = false, const bool only_explicit = false) const;
  void dump(const BNGData& bng_data, const bool for_diff, const bool only_explicit, const std::string ind = "") const;
};

typedef small_vector<MolInstance> MolInstanceVector;

} /* namespace BNG */

#endif /* LIBS_BNG_MOL_INSTANCE_H_ */
