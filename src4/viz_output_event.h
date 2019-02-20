/*
 * viz_output_event.h
 *
 *  Created on: Feb 19, 2019
 *      Author: adam
 */

#ifndef SRC4_VIZ_OUTPUT_EVENT_H_
#define SRC4_VIZ_OUTPUT_EVENT_H_

#include "base_event.h"

namespace mcell {

class viz_output_event_t: public base_event_t {
public:
	viz_output_event_t(world_t* world_)
		: base_event_t(EVENT_TYPE_INDEX_VIZ_OUTPUT),
			viz_mode(NO_VIZ_MODE), file_prefix_name(nullptr),
			world(world_) {
	}
	virtual ~viz_output_event_t() {}

	virtual void step();
	virtual void dump(const std::string indent);

	viz_mode_t viz_mode;
	const char* file_prefix_name; // in const pool

	world_t* world;
private:
	FILE* create_and_open_output_file_name();
	void output_ascii_molecules();
	void output_cellblender_molecules();
};

} /* namespace mcell */

#endif /* SRC4_VIZ_OUTPUT_EVENT_H_ */
