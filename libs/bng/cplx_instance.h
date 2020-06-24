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
#include "bng/graph.h"

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
  CplxInstance(const BNGData* bng_data_)
    : orientation(ORIENTATION_NONE),
      bng_data(bng_data_)
      {
  }

  CplxInstance(const CplxInstance& other)
    : mol_instances(other.mol_instances),
      orientation(other.orientation),
      bng_data(other.bng_data)
      {
    // copy ctor is needed because we must recreate graph that has pointers to
    // molecule and complex instances, finalize also sets flags
    finalize();
  }


  // must be called after initialization, sets up flags
  // also creates graphs for non-simple complexes
  void finalize();

  void create_graph();

  const Graph& get_graph() const {
    assert(is_finalized());
    return graph;
  }

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
    assert(is_finalized() && pattern.is_finalized());
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
#ifdef DEBUG_CPLX_MATCHING_EXTRA_COMPARE
      std::cout << "Comparing: ";
      dump(false);
      std::cout << " with: ";
      other.dump(false); std::cout << "\n";
#endif
    assert(is_finalized() && other.is_finalized());
    if (is_simple() && other.is_simple()) {
      return matches_simple(other, ignore_orientation);
    }
    else if (is_simple() != other.is_simple()) {
      // cannot be identical when
      return false;
    }
    return (!ignore_orientation || (orientation == other.orientation)) &&
        matches_complex_fully_ignore_orientation(other);
  }

  // full match & all other members must me equal
  bool operator == (const CplxInstance& other) const {
    assert(is_finalized() && other.is_finalized());

    // ordering of components in a molecule is not important
    return
        orientation == other.orientation &&
        get_flags() == other.get_flags() &&
        matches_fully(other);
  }

  std::string to_str(const BNGData& bng_data, bool in_reaction = false) const;
  void dump(const bool for_diff = false, const std::string ind = "") const;

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

  // boost's graph isomorphism does not cope with constant graphs
  // TODO: not sure of the internal representation really changes, might have
  // impact on parallel execution
  mutable Graph graph;

  const BNGData* bng_data; // needed mainly for dumps and debugging
};

typedef small_vector<CplxInstance> CplxInstanceVector;

} /* namespace BNG */

#endif /* LIBS_BNG_CPLX_INSTANCE_H_ */
