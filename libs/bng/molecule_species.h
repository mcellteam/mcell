/*
 * elementary_molecule.h
 *
 *  Created on: Jan 9, 2020
 *      Author: ahusar
 */

#ifndef LIBS_BNG_MOLECULE_SPECIES_H_
#define LIBS_BNG_MOLECULE_SPECIES_H_

#include <string>

#include "bng_defines.h"
#include "component.h"

namespace BNG {

// elementary molecules exist only once per instance
// a change in state of its component means that
// we have a different elementary molecule ID
class MoleculeSpecies {
public:
  std::string name;

  small_vector<Component> components;
};

} /* namespace BNG */

#endif /* LIBS_BNG_MOLECULE_SPECIES_H_ */
