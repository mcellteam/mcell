/*
 * elementary_molecule.h
 *
 *  Created on: Jan 9, 2020
 *      Author: ahusar
 */

#ifndef LIBS_BNG_ELEM_MOL_TYPE_H_
#define LIBS_BNG_ELEM_MOL_TYPE_H_

#include <string>

#include "bng/bng_defines.h"
#include "bng/base_flag.h"

namespace BNG {

class BNGData;

class ComponentType {
public:
  // name of the component itself is not sufficient to uniquely identify it,
  // it must be always searched with its elementary molecule type name in context
  std::string name;
  // using elementary molecule type name instead of ID becaue during construction
  std::string elem_mol_type_name;

  uint_set<state_id_t> allowed_state_ids;

  bool operator ==(const ComponentType& ct2) const {
    // ordering of allowed states is not important (they are in a set anyway)
    // two states must have the same id, this is ensured in find_or_add_state_name
    return name == ct2.name && allowed_state_ids == ct2.allowed_state_ids;
  }

  std::string to_str(const BNGData& bng_data) const;
  void dump(const BNGData& bng_data) const;
};


// Information common both to MolType and Species
class ElemMolTypeSpeciesCommonData {
public:
  ElemMolTypeSpeciesCommonData()
    : D(FLT_INVALID),
      time_step(FLT_INVALID), space_step(FLT_INVALID),
      color_set(false), color_r(1), color_g(0), color_b(0), scale(1) {
  }

  ElemMolTypeSpeciesCommonData(const ElemMolTypeSpeciesCommonData& other)
    : D(other.D),
      time_step(other.time_step),
      space_step(other.space_step),
      color_set(other.color_set),
      color_r(other.color_r), color_g(other.color_g), color_b(other.color_b),
      scale(other.scale) {
  }

  float_t D; // diffusion constant

  // computed values, set differently for ElemMolType and Species
  float_t time_step;
  float_t space_step;

  // visualization
  void set_color(float_t r, float_t g, float_t b) {
    color_r = r;
    color_g = g;
    color_b = b;
    color_set = true;
  }
  void set_scale(float_t s) {
    scale = s;
  }

  bool color_set;
  float_t color_r, color_g, color_b ;  // mol color default is red
  float_t scale; // scale = 1 by default

};


// Molecule type determines all allowed components and states of these components.
// It is only used to check that reactions and instantiations (releases) follow the
// allowed components and states.
class ElemMolType: public BaseSpeciesCplxMolFlag, public ElemMolTypeSpeciesCommonData {
public:
  ElemMolType() :
    custom_time_step(0), custom_space_step(0) {
  }

  std::string name;
  std::vector<component_type_id_t> component_type_ids;

  // when the user supplied a custom step, the attribute
  // is set to non-zero value, max one of them can be set to a non-zero value
  float_t custom_time_step;
  float_t custom_space_step;

  bool has_custom_time_or_space_step() const {
    return custom_time_step != 0 || custom_space_step != 0;
  }

  void compute_space_and_time_step(const BNGConfig& config);

  bool cant_initiate() const {
    return has_flag(SPECIES_MOL_FLAG_TARGET_ONLY);
  }

  bool operator ==(const ElemMolType& mt2) const {
    // ordering of components in a molecule is important
    // two component types must have the same id, this is ensured in find_or_add_component_type
    // diffusion constant and custom time/space steps are ignored,
    // flags are ignored as well
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

  std::string to_str(const BNGData& bng_data) const;
  void dump(const BNGData& bng_data) const;
};


// auxiliary functions used by ElemMolType and Species
float_t get_default_space_step(const BNGConfig& config, const float_t D);

void get_space_and_time_step(
    const BNGConfig& config,
    const bool is_surf, const float_t D,
    const float_t custom_time_step, const float_t custom_space_step,
    float_t& time_step, float_t& space_step);

} /* namespace BNG */

#endif /* LIBS_BNG_ELEM_MOL_TYPE_H_ */
