/*
 * complex_species.h
 *
 *  Created on: Jan 9, 2020
 *      Author: ahusar
 */

#ifndef LIBS_BNG_COMPLEX_SPECIES_H_
#define LIBS_BNG_COMPLEX_SPECIES_H_

#include "bng_defines.h"
#include "molecule_type.h"

namespace BNG {


class ComponentInstance {
public:
  component_type_id_t component_type_id;
  state_id_t state;
  bond_value_t bond_value;
};

// This extension also captures bonds and specific
// not using virtual methods
// both pattern and instance
// Similarly as ComplexSpeciesInstance -
//  can be used as instance and as pattern
class MoleculeTypeInstance {
public:
  // ID of this molecule type in BNGData::molecule_types
  molecule_type_id_t molecule_type_id;

  small_vector<ComponentInstance> component_instances;
};

// this class is used in two ways:
// - as a pattern for matching, not all states and bonds need to be entered
// - as a definition of species, in this case all components must be present and
//      if a component has more than 0 states then the state must be set
class ComplexSpeciesInstance {
public:
  small_vector<MoleculeTypeInstance> molecule_patterns;

  bool operator ==(const ComplexSpeciesInstance& cs2) const  {
    // ordering of components in a molecule is important
    // two component types must have the same id, this is ensured in find_or_add_component_type
    return molecule_patterns == cs2.molecule_patterns;
  }
};


// maybe some derived class for instances?

typedef small_vector<ComplexSpeciesInstance> ComplexSpeciesInstanceVector;

} /* namespace BNG */

#endif /* LIBS_BNG_COMPLEX_SPECIES_H_ */
