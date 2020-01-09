/*
 * component.h
 *
 *  Created on: Jan 9, 2020
 *      Author: ahusar
 */

#ifndef LIBS_BNG_COMPONENT_H_
#define LIBS_BNG_COMPONENT_H_

#include <string>

#include "bng_defines.h"

namespace BNG {

class ComponentType {
  std::string name;
  uint_set<state_id_t> allowed_states;
};


// each Component belongs directly to its elementary molecule
class Component {
public:
  component_type_id_t component_type_id;
  state_id_t state; // state must be one of the states in ComponentType::allowed_states

  // information on position of this component (binding site) will be here

  // info on binding is stored in ComplexSpecies
};

} /* namespace BNG */

#endif /* LIBS_BNG_COMPONENT_H_ */
