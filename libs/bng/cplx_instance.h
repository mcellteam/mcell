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
#include "base_flag.h"
#include "mol_instance.h"

// rename this to complex instance?

namespace BNG {

class MolType;
class BNGData;


// this class is used in two ways:
// - as a pattern for matching, not all states and bonds need to be entered
// - as a definition of species, in this case all components must be present and
//      if a component has more than 0 states then the state must be set
// TODO:L rename to pattern
class CplxInstance: public BaseFlag {
public:
  small_vector<MolInstance> mol_patterns;

private:
  // ID of this pattern, set to a value only
  // if this instance is used by a reaction
  // ???
  // pattern_id_t pattern_id;

  // not read from BNG yet, but proposal is on its way
  // for now
  orientation_t orientation;

public:
  CplxInstance()
    : orientation(ORIENTATION_NONE) {
  }

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

  /*bool has_single_orientation() const {
    return has_flag(CPLX_FLAG_HAS_SINGLE_ORIENTATION);
  }

  // asserts if has_single_orientation is false
  orientation_t get_single_orientation() const {
    assert(has_single_orientation());
    return has_flag(CPLX_FLAG_SINGLE_ORIENTATION_IS_UP) ? ORIENTATION_UP : ORIENTATION_DOWN;
  }*/

  bool is_simple() const {
    return has_flag(CPLX_FLAG_ONE_MOL_NO_COMPONENTS);
  }

  // asserts if has_single_orientation is false
  //void set_single_orientation(orientation_t orientation) const;

  bool get_orientation() const {
    return orientation;
  }

  void set_orientation(const orientation_t o) {
    orientation = o;
  }

  // returns true if this object as a pattern matches second instance
  bool matches(const CplxInstance& inst, const bool ignore_orientation = false) const;

  bool operator ==(const CplxInstance& ci2) const {
    // ordering of components in a molecule is important
    // two component types must have the same id, this is ensured in find_or_add_component_type
    return
        mol_patterns == ci2.mol_patterns &&
        orientation == ci2.orientation;
  }

  void dump(const BNGData& bng_data) const;
};


// maybe some derived class for instances?

typedef small_vector<CplxInstance> CplxInstanceVector;

} /* namespace BNG */

#endif /* LIBS_BNG_CPLX_INSTANCE_H_ */
