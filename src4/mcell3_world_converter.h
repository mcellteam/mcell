/*
 * mcell3_world_converter.h
 *
 *  Created on: Jan 30, 2019
 *      Author: adam
 */

#ifndef SRC4_MCELL3_WORLD_CONVERTER_H_
#define SRC4_MCELL3_WORLD_CONVERTER_H_

#include "mcell_structs.h"

#ifdef __cplusplus
extern "C"
#endif
bool mcell4_convert_mcell3_volume(struct volume* s);

#ifdef __cplusplus
extern "C"
#endif
bool mcell4_run_simulation();

#ifdef __cplusplus
extern "C"
#endif
void mcell4_delete_world();


// following code is only for C++
#ifdef __cplusplus

#include "world.h"
#include <map>

namespace mcell {

class mcell3_world_converter {
public:
	mcell3_world_converter() :
		world(nullptr) {
	}

	~mcell3_world_converter() {
		delete world;
	}

	void reset();

	bool convert(volume* s);
	bool convert_simulation_setup(volume* s);
	bool convert_species(volume* s);
	bool convert_release_events(volume* s);

	// contained world object
	world_t* world;

private:
	species_id_t get_new_species_id(u_int mcell3_id) {
		auto it = mcell3_species_id_map.find(mcell3_id);
		assert(it != mcell3_species_id_map.end());
		return it->second;
	}

	std::map<u_int, species_id_t> mcell3_species_id_map;
};


} // namespace mcell

#endif // #ifdef __cplusplus

#endif /* SRC4_MCELL3_WORLD_CONVERTER_H_ */
