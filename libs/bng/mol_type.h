/*
 * elementary_molecule.h
 *
 *  Created on: Jan 9, 2020
 *      Author: ahusar
 */

#ifndef LIBS_BNG_MOL_TYPE_H_
#define LIBS_BNG_MOL_TYPE_H_

#include <string>

#include "bng/bng_defines.h"
#include "bng/base_flag.h"

namespace BNG {

class BNGData;

class ComponentType {
public:
  std::string name;
  uint_set<state_id_t> allowed_state_ids;

  bool operator ==(const ComponentType& ct2) const {
    // ordering of allowed states is not important (they are in a set anyway)
    // two states must have the same id, this is ensured in find_or_add_state_name
    return name == ct2.name && allowed_state_ids == ct2.allowed_state_ids;
  }

  void dump(const BNGData& bng_data) const;
};


// Molecule type determines all allowed components and states of these components.
// It is only used to check that reactions and instantiations (releases) follow the
// allowed components and states.
class MolType: public BaseSpeciesCplxMolFlag {
public:
  MolType()
    : D(FLT_INVALID), custom_time_step(0), custom_space_step(0) {
  }

  std::string name;
  small_vector<component_type_id_t> component_type_ids;

  float_t D; // diffusion constant

  // when the user supplied a custom step, the attribute
  // is set to non-zero value, max one of them can be set to a non-zero value
  float_t custom_time_step;
  float_t custom_space_step;

  bool operator ==(const MolType& mt2) const {
    // ordering of components in a molecule is important
    // two component types must have the same id, this is ensured in find_or_add_component_type
    // diffusion constant and custom time/space steps are ignored
    return name == mt2.name && component_type_ids == mt2.component_type_ids;
  }

  uint get_component_uses_count(const component_type_id_t counted_ct_id) const {
    uint res = 0;
    for (component_type_id_t ct_id: component_type_ids) {
      if (counted_ct_id == ct_id) {
        res++;
      }
    }
    return res;
  }

  void dump(const BNGData& bng_data) const;
};

} /* namespace BNG */

#endif /* LIBS_BNG_MOL_TYPE_H_ */
