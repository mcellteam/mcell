/*
 * species.h
 *
 *  Created on: Feb 10, 2019
 *      Author: adam
 */

#ifndef SRC4_SPECIES_H_
#define SRC4_SPECIES_H_

#include <string>
#include "defines.h"

namespace mcell {

class species {
	uint16_t species_id;
	float_t D; // diffusion constant
	std::string name;
};

} /* namespace mcell */

#endif /* SRC4_SPECIES_H_ */
