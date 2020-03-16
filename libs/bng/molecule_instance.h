/*
 * base_molecule.h
 *
 *  Created on: Jan 9, 2020
 *      Author: ahusar
 */

#ifndef LIBS_BNG_MOLECULE_INSTANCE_H_
#define LIBS_BNG_MOLECULE_INSTANCE_H_

#include "bng_defines.h"

namespace BNG {


enum class BaseMoleculeFlags {

  // this is the single molecule of a complex that is being
  // actively diffused,
  // other molecules of this complex are used for detection of reactions,
  // but the diffusion of this specific molecule moves them all
  //
  // should have index 0 in the ComplexSpecies molecules array
  MainMolecule = 1 << 0,

  // this molecule belongs to a complex that is not represented by a
  // single point
  // when not set, reactions are detected in a way that the whole complexes
  // react against each other
  SpatialComplex = 1 << 1,


  // Temporary marker, equivalent to MCell 3 EXTERNAL_SPECIES
  // Will be removed in the future but for now we would like to support both MCell3
  // and BNG reactions to produce the same result as in MCell-R
  BNGSpecies = 1 << 31,
};

// base class for molecules,
// no virtual methods are allowed because molecules are essential
// objects that must be fast
class MoleculeInstance {
public:
  // to what complex species this molecule belongs
  complex_instance_index_t complex_instance_index;

  // index to the complex species's elementary molecules vector
  // that represents this specific instance
  uint elementary_molecule_index;

  // flags formed from BaseMoleculeFlags
  uint base_flags;
};

// holds information about a complex instance

// TODO: who owns and destroys complex instances? - BNGEngine
class ComplexInstance {
public:
  // to what complex species this complex belongs
  complex_species_id_t complex_species_index;

  // molecule instances contained in this complex instance
  small_vector<molecule_id_t> molecule_ids;

};

} /* namespace BNG */

#endif /* LIBS_BNG_MOLECULE_INSTANCE_H_ */
