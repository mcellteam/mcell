/*
 * complex_species.h
 *
 *  Created on: Jan 9, 2020
 *      Author: ahusar
 */

#ifndef LIBS_BNG_COMPLEX_SPECIES_H_
#define LIBS_BNG_COMPLEX_SPECIES_H_

#include "bng_defines.h"

// rename this to complex instance?

namespace BNG {

class MoleculeType;
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
  ComponentInstance(const component_type_id_t id)
    : component_type_id(id), state_id(STATE_ID_DONT_CARE), bond_value(BOND_VALUE_NO_BOND) {
  }

  // type id can be also found from parent's MoleculeInstance::molecule_type_id
  component_type_id_t component_type_id;

  // specific state or STATE_ID_DONT_CARE if we don't care
  state_id_t state_id;

  // number for bonds (0, 1, ...),
  // BOND_VALUE_ANY for patterns (+),
  // BOND_VALUE_NO_BOND if this component has no bond
  bond_value_t bond_value;

  bool operator ==(const ComponentInstance& comp2) const  {
    return
        component_type_id == comp2.component_type_id &&
        state_id == comp2.state_id &&
        bond_value == comp2.bond_value;
  }

  void dump(const BNGData& bng_data) const;
};


// Similarly as ComplexSpeciesInstance -
//  can be used as instance and as pattern
// TODO: should molecule type reference be an attribute?
class MoleculeInstance {
public:
  // ID of this molecule type in BNGData::molecule_types
  molecule_type_id_t molecule_type_id;

  // has the same number of elements as MoleculeType::component_type_ids
  small_vector<ComponentInstance> component_instances;


  bool operator ==(const MoleculeInstance& mi2) const  {
    return molecule_type_id == mi2.molecule_type_id && component_instances == mi2.component_instances;
  }


  void initialize_components_types(const MoleculeType& mt);

  // searches for component with name
  uint get_corresponding_component_index(
      const BNGData& bng_data,
      const MoleculeType& mt,
      const std::string& name,
      const uint starting_index
  ) const;

  void dump(const BNGData& bng_data) const;
};

// this class is used in two ways:
// - as a pattern for matching, not all states and bonds need to be entered
// - as a definition of species, in this case all components must be present and
//      if a component has more than 0 states then the state must be set
class ComplexInstance {
public:
  small_vector<MoleculeInstance> molecule_patterns;

  bool operator ==(const ComplexInstance& ci2) const  {
    // ordering of components in a molecule is important
    // two component types must have the same id, this is ensured in find_or_add_component_type
    return molecule_patterns == ci2.molecule_patterns;
  }

  void dump(const BNGData& bng_data) const;
};


// maybe some derived class for instances?

typedef small_vector<ComplexInstance> ComplexInstanceVector;

} /* namespace BNG */

#endif /* LIBS_BNG_COMPLEX_SPECIES_H_ */
