/*
 * complex_species.h
 *
 *  Created on: Jan 9, 2020
 *      Author: ahusar
 */

#ifndef LIBS_BNG_CPLX_H_
#define LIBS_BNG_CPLX_H_

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
class Cplx: public BaseSpeciesCplxMolFlag {
public:
  MolInstanceVector mol_instances;

private:
  // not read from BNG yet, but proposal is on its way
  // used in reactions, not in species
  orientation_t orientation;

  compartment_id_t compartment_id;

  // all possible compartments where this cplx might be used in reactions
  CompartmentIdSet reactant_compartments;

public:
  Cplx(const BNGData* bng_data_)
    : orientation(ORIENTATION_NONE),
      compartment_id(COMPARTMENT_ID_INVALID),
      bng_data(bng_data_)
      {
  }

  Cplx(const Cplx& other) {
    *this = other;
  }

  Cplx& operator =(const Cplx& other) {
    mol_instances = other.mol_instances;
    orientation = other.orientation;
    compartment_id = other.compartment_id;
    bng_data = other.bng_data;

    set_flags(other.get_flags());

    // copy ctor is needed because we must recreate graph that has pointers to
    // molecule and complex instances, finalize also sets some flags
    finalize();

    return *this;
  }


  // must be called after initialization, sets up flags
  // also creates graphs for non-simple complexes
  void finalize();

  // called automatically from finalize but needed elsewhere as well
  void update_flag_and_compartments_used_in_rxns();

  void create_graph();

  const Graph& get_graph() const {
    assert(is_finalized());
    return graph;
  }

  Graph& get_graph() {
    assert(is_finalized());
    return graph;
  }

  bool is_simple() const {
    return has_flag(SPECIES_CPLX_FLAG_ONE_MOL_NO_COMPONENTS);
  }

  bool is_canonical() const {
    return has_flag(SPECIES_CPLX_FLAG_IS_CANONICAL);
  }

  // returns true if all components of all molecules are be present and their is state set
  bool is_fully_qualified() const;

  // returns true if all molecules are connected through bonds,
  // this is the only currently supported mode and it is checked in semantics analyzer
  bool is_connected() const;

  mol_type_id_t get_simple_species_mol_type_id() const {
    assert(is_simple());
    return mol_instances[0].mol_type_id;
  }

  orientation_t get_orientation() const {
    return orientation;
  }

  void set_orientation(const orientation_t o) {
    // TODO: here could be some extra checks related to compartments
    orientation = o;
  }

  bool has_compartment() const {
    assert(compartment_id != COMPARTMENT_ID_INVALID);
    return compartment_id != COMPARTMENT_ID_NONE && compartment_id != COMPARTMENT_ID_ANY;
  }

  compartment_id_t get_compartment_id() const {
    return compartment_id;
  }

  void set_compartment_id(const compartment_id_t cid) {
    // TODO: here could be some extra checks related to orientation
    compartment_id = cid;
  }

  // returns true if this object as a pattern matches second instance
  bool matches_pattern(const Cplx& pattern, const bool ignore_orientation = false) const {
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

  // returns how many times a pattern matches this cplx instance,
   // used for counting of molecule pattern type observables
  uint get_pattern_num_matches(const Cplx& pattern) const;

  // returns true if this complex is equivalent
  bool matches_fully(const Cplx& other, const bool ignore_orientation = false) const {
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
      // cannot be identical when one is simple and the other not
      return false;
    }
    return (!ignore_orientation || (orientation == other.orientation)) &&
        matches_complex_fully_ignore_orientation(other);
  }

  // full match & all other members must me equal
  bool operator == (const Cplx& other) const {
    assert(is_finalized() && other.is_finalized());

    // ordering of components in a molecule is not important
    return
        orientation == other.orientation &&
        get_flags() == other.get_flags() &&
        matches_fully(other);
  }

  // sort molecule instances so that all identical complexes use
  // the same ordering
  void canonicalize();

  const CompartmentIdSet& get_reactant_compartments() const {
    assert(!reactant_compartments.empty() && "Not initialized");
    return reactant_compartments;
  }

  bool is_reactant_compartment(const compartment_id_t compartment_id) const {
    assert(!reactant_compartments.empty() && "Not initialized");
    return reactant_compartments.count(compartment_id) != 0;
  }

  // if the compartment passed as argument is in a set of applicable compartments then
  // returns its value, otherwise returns COMPARTMENT_ID_NONE because we must not
  // set a compartment ID that is not applicable, e.g. the RxnContaines counts on it
  compartment_id_t get_as_reactant_compartment(const compartment_id_t compartment_id) const {
    assert(!reactant_compartments.empty() && "Not initialized");
    if (is_reactant_compartment(compartment_id)) {
      return compartment_id;
    }
    else {
      // outside of any compartment important for this species
      return COMPARTMENT_ID_NONE;
    }
  }

  std::string to_str(bool in_surf_reaction = false) const;
  void dump(const bool for_diff = false, const std::string ind = "") const;

private:
  bool matches_simple(const Cplx& other, const bool ignore_orientation = false) const {
    assert(is_simple() && other.is_simple());
    assert(mol_instances.size() == 1 && other.mol_instances.size() == 1);

    if (!ignore_orientation && orientation != other.orientation) {
      return false;
    }
    return mol_instances[0].matches_simple(other.mol_instances[0]);
  }

  bool matches_complex_pattern_ignore_orientation(const Cplx& pattern) const;
  bool matches_complex_fully_ignore_orientation(const Cplx& pattern) const;

  void sort_components_and_mols();
  void renumber_bonds();

  // boost's graph isomorphism does not cope with constant graphs
  // TODO: not sure of the internal representation really changes, might have
  // impact on parallel execution
  mutable Graph graph;

  const BNGData* bng_data; // needed mainly for dumps and debugging
};

typedef small_vector<Cplx> CplxVector;

} /* namespace BNG */

#endif /* LIBS_BNG_CPLX_H_ */
