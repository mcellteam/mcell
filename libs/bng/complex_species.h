/*
 * complex_species.h
 *
 *  Created on: Jan 9, 2020
 *      Author: ahusar
 */

#ifndef LIBS_BNG_COMPLEX_SPECIES_H_
#define LIBS_BNG_COMPLEX_SPECIES_H_

#include "bng_defines.h"
#include "molecule_type.h"

namespace BNG {

// a bond is simply a pair of indices to ComplexSpecies::molecule_patterns
//   (or equivalent ComplexSpecies::molecule_ids)
// represents an edge in the molecule graph
class Bond {
public:
  uint index1;
  uint index2;

  bool operator ==(const Bond& b2) const {
    return index1 == b2.index1 && index2 == b2.index2;
  }
};


// this class is used in two ways:
// - as a pattern for matching, not all states and bonds need to be entered
// - as a definition of species, in this case all components must be present and
//      if a component has more than 0 states then the state must be set
//
// but still, we need???
class ComplexSpecies {
public:
  // the two arrays must have the same size, molecule_ids serve for fast comparison
  small_vector<molecule_type_id_t> molecule_ids;
  small_vector<MoleculeType> molecule_patterns;

  // bonds - must be defined in a way that allows for fast graph
  // algorithms
  small_vector<Bond> bonds;

  //float_t diffusion_constant;
  bool operator ==(const ComplexSpecies& cs2) const  {
    assert(molecule_ids.size() == molecule_patterns.size());
    // ordering of components in a molecule is important
    // two component types must have the same id, this is ensured in find_or_add_component_type
    return molecule_ids == cs2.molecule_ids && molecule_patterns == cs2.molecule_patterns && bonds == cs2.bonds;
  }
};


// maybe some derived class for instances?

typedef small_vector<ComplexSpecies> ComplexSpeciesVector;

} /* namespace BNG */

#endif /* LIBS_BNG_COMPLEX_SPECIES_H_ */
