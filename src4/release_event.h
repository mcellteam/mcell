/*
 * release_molecules.h
 *
 *  Created on: Jan 30, 2019
 *      Author: adam
 */

#ifndef SRC4_RELEASE_EVENT_H_
#define SRC4_RELEASE_EVENT_H_

#include "base_event.h"

namespace mcell {

class release_event_t: public base_event_t {
public:
	release_event_t() :
		species_id(SPECIES_ID_INVALID),
		release_number(0)	{
	}
	virtual ~release_event_t() {}

	virtual void step() {}
	virtual void dump(const std::string ind);

	vec3_t location;
	species_id_t species_id;
	uint32_t release_number; // number of molecules to release
	std::string name;
};

} // namespace mcell


#endif /* SRC4_RELEASE_EVENT_H_ */
