/*
 * complex_species.h
 *
 *  Created on: Jan 9, 2020
 *      Author: ahusar
 */

#ifndef LIBS_BNG_CPLX_INSTANCE_H_
#define LIBS_BNG_CPLX_INSTANCE_H_

#include <iostream>

#include "bng_defines.h"

// rename this to complex instance?

namespace BNG {

class MolType;
class BNGData;


// NOTE: maybe use bitfield instead?
class BaseFlag {
private:
  bool finalized;
  uint flags;

public:
  BaseFlag()
    : finalized(false), flags(0) {
  }

  bool has_flag(uint flag) const {
    assert(finalized);
    assert(__builtin_popcount(flag) == 1);
    return (flags & flag) != 0;
  }

  void set_flag(uint flag, bool value = true) {
    assert(finalized);
    assert(__builtin_popcount(flag) == 1);
    if (value) {
      flags = flags | flag;
    }
    else {
      clear_flag(flag);
    }
  }

  void clear_flag(uint flag) {
    assert(finalized);
    assert(__builtin_popcount(flag) == 1);
    flags = flags & ~flag;
  }

  void set_finalized() {
    finalized = true;
  }
};

/**
 * Pattern vs instance:
 * - patterns are used in reaction rules, they are written in the don't care don't write manner and
 *   may use '+' as bond
 *
 * - instances are fully specified, each component has a given state (if it has states)
 */

class ComponentInstance {
public:
  ComponentInstance(const component_type_id_t id)
    : component_type_id(id), state_id(STATE_ID_DONT_CARE), bond_value(BOND_VALUE_ANY),
      explicitly_listed_in_pattern(false) {
  }

  // type id can be also found from parent's MoleculeInstance::molecule_type_id
  // maybe remove
  component_type_id_t component_type_id;

  // specific state or STATE_ID_DONT_CARE if we don't care
  state_id_t state_id;

  // number for bonds (0, 1, ...),
  // BOND_VALUE_ANY for patterns (+),
  // BOND_VALUE_NO_BOND if this component has no bond
  bond_value_t bond_value;

  // true if this component was explicitly listed in a pattern
  bool explicitly_listed_in_pattern;

  bool operator ==(const ComponentInstance& comp2) const  {
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

  void dump(const BNGData& bng_data, std::ostream& out = std::cout) const;
};


enum cplx_mol_flag_t {
  CPLX_MOL_FLAG_SURF = 1 << 0,
  CPLX_MOL_FLAG_VOL = 1 << 1,

  CPLX_FLAG_HAS_SINGLE_ORIENTATION = 1 << 2,
  CPLX_FLAG_SINGLE_ORIENTATION_IS_UP = 1 << 3,
  CPLX_FLAG_ONE_MOL_NO_COMPONENTS = 1 << 4,
};

// Similarly as ComplexSpeciesInstance -
//  can be used as instance and as pattern
// TODO: should molecule type reference be an attribute?
// TODO: rename to pattern
class MolInstance: public BaseFlag {
public:
  MolInstance()
    : mol_type_id(MOL_TYPE_ID_INVALID), orientation(ORIENTATION_NONE) {
  }

  // share this interface with actual species?
  bool is_vol() const {
    return has_flag(CPLX_MOL_FLAG_VOL);
  }

  // if any of the contained molecule instances is a surface molecule,
  // the whole complex is a surface molecule
  bool is_surf() const {
    return has_flag(CPLX_MOL_FLAG_SURF);
  }

  void finalize() {
    set_finalized();
    // flag about molecule type must be set
    assert(is_vol() || is_surf());
  }


  // ID of this molecule type in BNGData::molecule_types
  mol_type_id_t mol_type_id;

  // has the same number of elements as MoleculeType::component_type_ids
  small_vector<ComponentInstance> component_instances;

  // not read from BNG yet, but
  orientation_t orientation;

  bool operator ==(const MolInstance& mi2) const  {
    return
        mol_type_id == mi2.mol_type_id &&
        component_instances == mi2.component_instances &&
        orientation == mi2.orientation;
  }


  void initialize_components_types(const MolType& mt);

  // searches for component with name
  uint get_corresponding_component_index(
      const BNGData& bng_data,
      const MolType& mt,
      const std::string& name,
      const uint starting_index
  ) const;

  void dump(const BNGData& bng_data, const bool only_explicit = false, std::ostream& out = std::cout) const;
};


// this class is used in two ways:
// - as a pattern for matching, not all states and bonds need to be entered
// - as a definition of species, in this case all components must be present and
//      if a component has more than 0 states then the state must be set
// TODO:L rename to pattern
class CplxInstance: public BaseFlag {
public:
  small_vector<MolInstance> mol_patterns;

  // must be called after initialization, sets up flags
  void finalize();

  // share this interface with actual species?
  bool is_vol() const {
    return has_flag(CPLX_MOL_FLAG_VOL);
  }

  // if any of the contained molecule instances is a surface molecule,
  // the whole complex is a surface molecule
  bool is_surf() const {
    return has_flag(CPLX_MOL_FLAG_SURF);
  }

  bool has_single_orientation() const {
    return has_flag(CPLX_FLAG_HAS_SINGLE_ORIENTATION);
  }

  // asserts if has_single_orientation is false
  orientation_t get_single_orientation() const {
    assert(has_single_orientation());
    return has_flag(CPLX_FLAG_SINGLE_ORIENTATION_IS_UP) ? ORIENTATION_UP : ORIENTATION_DOWN;
  }

  bool is_simple() const {
    return has_flag(CPLX_FLAG_ONE_MOL_NO_COMPONENTS);
  }

  // asserts if has_single_orientation is false
  void set_single_orientation(orientation_t orientation) const;

  bool operator ==(const CplxInstance& ci2) const {
    // ordering of components in a molecule is important
    // two component types must have the same id, this is ensured in find_or_add_component_type
    return mol_patterns == ci2.mol_patterns;
  }

  void dump(const BNGData& bng_data) const;
};


// maybe some derived class for instances?

typedef small_vector<CplxInstance> CplxInstanceVector;

} /* namespace BNG */

#endif /* LIBS_BNG_CPLX_INSTANCE_H_ */
