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

class base_molecule_t {
public:
	base_molecule_t();

	uint16_t flags; // defunct by default, use bitfield instead?
	uint16_t species_index;
};


class volume_molecule_t : public base_molecule_t {
	vec3_t pos;
};

} // namespace mcell

#endif /* SRC4_MOLECULE_H_ */
