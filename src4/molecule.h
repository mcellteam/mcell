/*
 * molecule.h
 *
 *  Created on: Feb 10, 2019
 *      Author: adam
 */

#ifndef SRC4_MOLECULE_H_
#define SRC4_MOLECULE_H_

#include "defines.h"

namespace mcell {

class base_molecule {
public:
	base_molecule();

	uint16_t flags; // defunct by default, use bitfield instead?
	uint16_t species_index;
};


class volume_molecule : public base_molecule {
	vec3_t pos;
};

} // namespace mcell

#endif /* SRC4_MOLECULE_H_ */
