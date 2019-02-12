/*
 * mcell3_world_converter.h
 *
 *  Created on: Jan 30, 2019
 *      Author: adam
 */

#ifndef SRC4_MCELL3_WORLD_CONVERTER_H_
#define SRC4_MCELL3_WORLD_CONVERTER_H_

#include "world.h"

struct volume; // in mcell_structs.h

#ifdef __cplusplus
extern "C"
#endif
bool convert_mcell3_volume_to_mcell3_world(volume* s);


namespace mcell {

class mcell3_world_converter {
public:
	mcell3_world_converter() :
		world(nullptr) {
	}

	~mcell3_world_converter() {
		delete world;
	}



	bool convert(volume* s);
	bool convert_species(volume* s);
	bool convert_release_events(volume* s);

	// global object
	world_t* world;
};



} // namespace mcell


#endif /* SRC4_MCELL3_WORLD_CONVERTER_H_ */
