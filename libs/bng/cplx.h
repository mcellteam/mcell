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
#include "bng/graph.h"
#include "bng/elem_mol.h"

namespace BNG {

class ElemMolType;
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
  ElemMolVector elem_mols;

private:
  // not read from BNG yet, but proposal is on its way
  // used in reactions, not in species
  orientation_t orientation;

public:
  Cplx(const BNGData* bng_data_)
    : orientation(ORIENTATION_NONE),
      bng_data(bng_data_)
      {
  }

  Cplx(const Cplx& other) {
    *this = other;
  }

  Cplx& operator =(const Cplx& other) {
    elem_mols = other.elem_mols;
    orientation = other.orientation;
    bng_data = other.bng_data;

    set_flags(other.get_flags());

    // copy ctor is needed because we must recreate graph that has pointers to
    // molecule and complex instances, finalize also sets some flags
    if (other.is_finalized()) {
      finalize_cplx(false);
    }

    return *this;
  }


  // must be called after initialization, sets up flags
  // also creates graphs for non-simple complexes
  void finalize_cplx(const bool init_flags_and_compartments = true);

  void create_graph();

  // if dont_know_elem_mol_types is true, we do not know the types elementary molecules
  // (whether they are vol or surf),
  compartment_id_t get_primary_compartment_id(const bool dont_know_elem_mol_types = false) const {
    assert(is_finalized());
    if (elem_mols.size() == 1) {
      return elem_mols[0].compartment_id;
    }
    else {
      return get_complex_compartment_id(dont_know_elem_mol_types);
    }
  }

  void get_used_compartments(uint_set<compartment_id_t>& compartments) const;


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

  elem_mol_type_id_t get_simple_species_mol_type_id() const {
    assert(is_simple());
    return elem_mols[0].elem_mol_type_id;
  }

  orientation_t get_orientation() const {
    return orientation;
  }

  void set_orientation(const orientation_t o) {
    // TODO: here could be some extra checks related to compartments
    orientation = o;
  }

  bool has_compartment() const {
    compartment_id_t id = get_primary_compartment_id();
    assert(id != COMPARTMENT_ID_INVALID);
    return id != COMPARTMENT_ID_NONE;
  }

  bool has_compartment_class_in_out() const {
    compartment_id_t id = get_primary_compartment_id();
    assert(id != COMPARTMENT_ID_INVALID);
    return is_in_out_compartment_id(id);;
  }

  // sets compartment to all contained elementary molecules
  void set_compartment_id(const compartment_id_t cid, const bool override_only_compartment_none = false);

  // go trough all elementary molecules and if compartment is set to cid, set it to NONE
  void remove_compartment_from_elem_mols(const compartment_id_t cid);

  // returns true if this object as a pattern matches second instance
  bool matches_pattern(const Cplx& pattern, const bool ignore_orientation = false) const {
    assert(is_finalized() && pattern.is_finalized());
    if (is_simple() && pattern.is_simple()) {
      return matches_simple_pattern(pattern, ignore_orientation);
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

  // returns true if this complex is equivalent, neither of the complexes is a pattern
  bool matches_fully(const Cplx& other, const bool ignore_orientation = false) const {
#ifdef DEBUG_CPLX_MATCHING_EXTRA_COMPARE
      std::cout << "Comparing: ";
      dump(false);
      std::cout << " with: ";
      other.dump(false); std::cout << "\n";
#endif
    assert(is_finalized() && other.is_finalized());
    if (is_simple() && other.is_simple()) {
      return matches_simple_fully(other, ignore_orientation);
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
  // default sorting of components is according to molecule types
  void canonicalize(const bool sort_components_by_name_do_not_finalize = false);

  // appends to string res
  void to_str(std::string& res, const bool in_surf_reaction = false) const;

  std::string to_str(const bool in_surf_reaction = false) const;
  void dump(const bool for_diff = false, const std::string ind = "") const;

private:
  bool matches_simple_pattern(const Cplx& pattern, const bool ignore_orientation = false) const {
    assert(is_simple() && pattern.is_simple());
    assert(elem_mols.size() == 1 && pattern.elem_mols.size() == 1);

    if (!ignore_orientation && orientation != pattern.orientation) {
      return false;
    }
    return pattern.elem_mols[0].matches_simple_pattern(elem_mols[0]);
  }

  bool matches_simple_fully(const Cplx& pattern, const bool ignore_orientation = false) const {
    assert(is_simple() && pattern.is_simple());
    assert(elem_mols.size() == 1 && pattern.elem_mols.size() == 1);

    if (!ignore_orientation && orientation != pattern.orientation) {
      return false;
    }
    return pattern.elem_mols[0].matches_simple_fully(elem_mols[0]);
  }

  compartment_id_t get_complex_compartment_id(const bool override_is_surface_cplx = false) const;

  bool matches_complex_pattern_ignore_orientation(const Cplx& pattern) const;
  bool matches_complex_fully_ignore_orientation(const Cplx& other) const;

  void sort_components_and_mols();
  void renumber_bonds();

  // boost's graph isomorphism does not cope with constant graphs
  // TODO: not sure of the internal representation really changes, might have
  // impact on parallel execution
  mutable Graph graph;
protected:
  // needed for computation of time/space step in Species and for dumps and debugging
  const BNGData* bng_data;
};

typedef small_vector<Cplx> CplxVector;

} /* namespace BNG */

#endif /* LIBS_BNG_CPLX_H_ */
