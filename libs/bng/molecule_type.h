/*
 * elementary_molecule.h
 *
 *  Created on: Jan 9, 2020
 *      Author: ahusar
 */

#ifndef LIBS_BNG_MOLECULE_TYPE_H_
#define LIBS_BNG_MOLECULE_TYPE_H_

#include <string>

#include "bng_defines.h"
#include "component.h"

namespace BNG {

class ComponentType {
public:
  std::string name;
  uint_set<state_id_t> allowed_state_ids;

  bool operator ==(const ComponentType& ct2) {
    // ordering of allowed states is not important (they are in a set anyway)
    // two states must have the same id, this is ensured in find_or_add_state_name
    return name == ct2.name && allowed_state_ids == ct2.allowed_state_ids;
  }
};


// Molecule type determines all allowed components and states of these components.
// It is only used to check that reactions and instantiations (releases) follow the
// allowed components and states.
class MoleculeType {
public:
  std::string name;
  small_vector<component_type_id_t> component_type_ids;

  bool operator ==(const MoleculeType& mt2) {
    // ordering of components in a molecule is important
    // two component types must have the same id, this is ensured in find_or_add_component_type
    return name == mt2.name && component_type_ids == mt2.component_type_ids;
  }
};

} /* namespace BNG */

#endif /* LIBS_BNG_MOLECULE_TYPE_H_ */
