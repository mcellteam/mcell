/*
 * complex_species.h
 *
 *  Created on: Jan 9, 2020
 *      Author: ahusar
 */

#ifndef LIBS_BNG_CPLX_INSTANCE_H_
#define LIBS_BNG_CPLX_INSTANCE_H_

#include <iostream>

#include "bng/bng_defines.h"
#include "bng/base_flag.h"
#include "bng/mol_instance.h"

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
class CplxInstance: public BaseSpeciesCplxMolFlag {
public:
  MolInstanceVector mol_instances;

  // NOTE: Some higher-level bond information will be needed here

private:
  // not read from BNG yet, but proposal is on its way
  // used in reactions, not in species
  orientation_t orientation;

public:
  CplxInstance()
    : orientation(ORIENTATION_NONE) {
  }

  // must be called after initialization, sets up flags
  void finalize_flags();

  bool is_simple() const {
    return has_flag(SPECIES_CPLX_FLAG_ONE_MOL_NO_COMPONENTS);
  }

  mol_type_id_t get_simple_species_mol_type_id() const {
    assert(is_simple());
    return mol_instances[0].mol_type_id;
  }

  orientation_t get_orientation() const {
    return orientation;
  }

  void set_orientation(const orientation_t o) {
    orientation = o;
  }

  // returns true if this object as a pattern matches second instance
  bool matches_pattern(const CplxInstance& pattern, const bool ignore_orientation = false) const {
    if (is_simple() && pattern.is_simple()) {
      return matches_simple(pattern, ignore_orientation);
    }
    else {
      if (!ignore_orientation && orientation != pattern.orientation) {
        return false;
      }
      else {
        return matches_complex_pattern_ignore_orientation(pattern);
      }
    }
  }

  // returns true if this complex is equivalent
  bool matches_fully(const CplxInstance& other, const bool ignore_orientation = false) const {
    if (is_simple() && other.is_simple()) {
      return matches_simple(other, ignore_orientation);
    }
    return
        mol_instances == other.mol_instances && // XXX this must be a graph comparison
        (!ignore_orientation || (orientation == other.orientation));
  }

  // full match & all other members must me equal
  bool operator ==(const CplxInstance& other) const {
    // ordering of components in a molecule is not important

    return
        orientation == other.orientation &&
        get_flags() == other.get_flags() &&
        matches_fully(other);
  }

  std::string to_str(const BNGData& bng_data, bool in_reaction = false) const;
  void dump(const BNGData& bng_data, const bool for_diff, const std::string ind = "") const;

private:
  bool matches_simple(const CplxInstance& other, const bool ignore_orientation = false) const {
    assert(is_simple() && other.is_simple());
    assert(mol_instances.size() == 1 && other.mol_instances.size() == 1);

    if (!ignore_orientation && orientation != other.orientation) {
      return false;
    }
    return mol_instances[0].matches_simple(other.mol_instances[0]);
  }

  bool matches_complex_pattern_ignore_orientation(const CplxInstance& pattern) const;
  bool matches_complex_fully_ignore_orientation(const CplxInstance& pattern) const;

};

typedef small_vector<CplxInstance> CplxInstanceVector;

} /* namespace BNG */

#endif /* LIBS_BNG_CPLX_INSTANCE_H_ */
