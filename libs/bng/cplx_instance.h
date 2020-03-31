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

namespace BNG {

class MolType;
class BNGData;


/**
 * Complex instance or pattern.
 *
 * This class is used in two ways:
 * - as a pattern for matching, not all states and bonds need to be entered
 * - as a definition of species, in this case all components must be present and
 *      if a component has more than 0 states then the state must be set
 *
 * Contains information on orientation, so two identical complexes but
 * with different orientation are different species.
 */
class CplxInstance: public BaseFlag {
public:
  MolInstanceVector mol_instances;

private:
  // not read from BNG yet, but proposal is on its way
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

  bool is_simple() const {
    return has_flag(CPLX_FLAG_ONE_MOL_NO_COMPONENTS);
  }

  orientation_t get_orientation() const {
    return orientation;
  }

  void set_orientation(const orientation_t o) {
    orientation = o;
  }

  // returns true if this object as a pattern matches second instance
  bool matches(const CplxInstance& inst, const bool ignore_orientation = false) const;

  bool equal(const CplxInstance& ci2) const {
    return
        mol_instances == ci2.mol_instances &&
        orientation == ci2.orientation;
  }

  bool operator ==(const CplxInstance& ci2) const {
    // ordering of components in a molecule is important
    // two component types must have the same id, this is ensured in find_or_add_component_type
    return equal(ci2);
  }

  void dump(const BNGData& bng_data, std::string ind = "") const;
};

typedef small_vector<CplxInstance> CplxInstanceVector;

} /* namespace BNG */

#endif /* LIBS_BNG_CPLX_INSTANCE_H_ */
