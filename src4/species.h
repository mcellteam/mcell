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

class species_t {
public:
	species_id_t species_id;

	uint32_t mcell3_species_id;
	float_t D; // diffusion constant
	std::string name;
	float_t space_step;
	float_t time_step;

	void dump(const std::string ind);
};

} /* namespace mcell */

#endif /* SRC4_SPECIES_H_ */
