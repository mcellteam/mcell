/*
 * complex_species.h
 *
 *  Created on: Jan 9, 2020
 *      Author: ahusar
 */

#ifndef LIBS_BNG_COMPLEX_SPECIES_H_
#define LIBS_BNG_COMPLEX_SPECIES_H_

#include "bng_defines.h"
#include "molecule_species.h"

namespace BNG {

// a bond is simply a pair of indices to ComplexSpecies::molecules
// represents an undirectional edge in the molecule graph
class Bond {
public:
  uint index1;
  uint index2;
};


class ComplexSpecies {
public:
  complex_species_index_t index;

  small_vector<MoleculeSpecies> molecules;

  // bonds - must be defined in a way that allows for fast graph
  // algorithms
  small_vector<Bond> bonds;


  float_t diffusion_constant;
};

} /* namespace BNG */

#endif /* LIBS_BNG_COMPLEX_SPECIES_H_ */
