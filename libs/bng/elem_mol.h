/*
 * complex_species.h
 *
 *  Created on: Jan 9, 2020
 *      Author: ahusar
 */

#ifndef LIBS_BNG_ELEM_MOL_H_
#define LIBS_BNG_ELEM_MOL_H_

#include <iostream>

#include "bng/bng_defines.h"
#include "bng/base_flag.h"

namespace BNG {

class ElemMolType;
class BNGData;

/**
 * Pattern vs instance:
 * - patterns are used in reaction rules, they are written in the don't care don't write manner and
 *   may use '+' as bond
 *
 * - instances are fully specified, each component has a given state (if it has states)
 */
class Component {
public:
  Component() :
    component_type_id(COMPONENT_TYPE_ID_INVALID), state_id(STATE_ID_DONT_CARE), bond_value(BOND_VALUE_BOUND) {
  }

  Component(const component_type_id_t id) :
    component_type_id(id), state_id(STATE_ID_DONT_CARE), bond_value(BOND_VALUE_BOUND) {
  }

  // type id can be also found from parent's MoleculeInstance::molecule_type_id
  // maybe remove
  component_type_id_t component_type_id;

  // specific state or STATE_ID_DONT_CARE if we don't care or unset
  state_id_t state_id;

  // number for bonds (0, 1, ...),
  // BOND_VALUE_BOUND for patterns (+),
  // BOND_VALUE_ANY for patterns (?),
  // BOND_VALUE_UNBOND if this component has no bond or a
  //    when used in a pattern that it must not have a bond
  //    (in cases where we do not care, the component is not listed at all)
  bond_value_t bond_value;

  bool operator == (const Component& comp2) const  {
    return
        component_type_id == comp2.component_type_id &&
        state_id == comp2.state_id &&
        bond_value == comp2.bond_value;
  }

  bool bond_has_numeric_value() const {
    assert(bond_value != BOND_VALUE_INVALID);
    return
        bond_value != BOND_VALUE_UNBOUND &&
        bond_value != BOND_VALUE_BOUND &&
        bond_value != BOND_VALUE_ANY;
  }

  bool state_is_set() const {
    assert(state_id != STATE_ID_INVALID);
    return state_id != STATE_ID_DONT_CARE;
  }

  void to_str(const BNGData& bng_data, std::string& res) const;
  std::string to_str(const BNGData& bng_data) const;
  void dump(const BNGData& bng_data, const std::string& ind = "") const;
};


/**
 * Instance of an elementary molecule type.
 * Named molecule type in BNGL but this might be confusing because molecules are whole complexes.
 * Similarly as ComplexSpeciesInstance, it can be used as instance and as pattern.
 */
class ElemMol: public BaseSpeciesCplxMolFlag {
public:
  // ID of this molecule type in BNGData::molecule_types
  elem_mol_type_id_t elem_mol_type_id;

  compartment_id_t compartment_id;

  // has the same number of elements as MoleculeType::component_type_ids
  small_vector<Component> components;

public:
  ElemMol() :
    elem_mol_type_id(MOL_TYPE_ID_INVALID), compartment_id(COMPARTMENT_ID_NONE)  {
  }

  // returns true if all components are present and their is state set
  bool is_fully_qualified(const BNGData& bng_data) const;

  void canonicalize(const BNGData& bng_data);

  // sorting by name is useful when we have no fixed molecule type
  // to guide the sorting
  void sort_components_by_name(const BNGData& bng_data);

  void finalize_flags_and_sort_components(const BNGData& bng_data);

  void insert_missing_components_as_any_state_pattern(const BNGData& bng_data);

  bool has_bond() const;

  // returns true if this object as a pattern matches second instance,
  // order is important for compartment matching
  bool matches_simple_pattern(const ElemMol& inst) const {
    assert(components.size() == 0 && inst.components.size() == 0 &&
        "Method can be used only for simple complexes, i.e. without components.");

    return elem_mol_type_id == inst.elem_mol_type_id &&
        (compartment_id == COMPARTMENT_ID_NONE || is_in_out_compartment_id(compartment_id) || compartment_id == inst.compartment_id);
  }

  bool matches_simple_fully(const ElemMol& inst) const {
    assert(components.size() == 0 && inst.components.size() == 0 &&
        "Method can be used only for simple complexes, i.e. without components.");

    return elem_mol_type_id == inst.elem_mol_type_id && compartment_id == inst.compartment_id;
  }

  bool operator == (const ElemMol& other) const  {
    return
        elem_mol_type_id == other.elem_mol_type_id &&
        compartment_id == other.compartment_id &&
        components == other.components &&
        get_flags() == other.get_flags();
  }

  // appends to string res
  void to_str(const BNGData& bng_data, std::string& res, const bool include_compartment = true) const;

  std::string to_str(const BNGData& bng_data) const;
  void dump(const BNGData& bng_data, const bool for_diff, const std::string ind = "") const;
};

typedef small_vector<ElemMol> ElemMolVector;

} /* namespace BNG */

#endif /* LIBS_BNG_ELEM_MOL_H_ */
