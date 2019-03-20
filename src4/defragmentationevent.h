/*
 * defragmentationevent.h
 *
 *  Created on: Mar 20, 2019
 *      Author: adam
 */

#ifndef SRC4_DEFRAGMENTATIONEVENT_H_
#define SRC4_DEFRAGMENTATIONEVENT_H_

#include "base_event.h"

namespace mcell {

class defragmentation_event: public base_event_t {
	defragmentation_event(world_t* world_)
		: base_event_t(EVENT_TYPE_INDEX_VIZ_OUTPUT),
			world(world_) {
	}

	virtual ~defragmentation_event() {};
	virtual void step();
	virtual void dump(const std::string indent) {};

	world_t* world;
};

} /* namespace mcell */

#endif /* SRC4_DEFRAGMENTATIONEVENT_H_ */
