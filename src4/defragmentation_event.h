/*
 * defragmentationevent.h
 *
 *  Created on: Mar 20, 2019
 *      Author: adam
 */

#ifndef SRC4_DEFRAGMENTATION_EVENT_H_
#define SRC4_DEFRAGMENTATION_EVENT_H_

#include "base_event.h"

namespace mcell {

class defragmentation_event_t: public base_event_t {
public:
	defragmentation_event_t(world_t* world_)
		: base_event_t(EVENT_TYPE_INDEX_DEFRAGMENTATION),
			world(world_) {
	}

	virtual void step();
	virtual void dump(const std::string indent);
private:
	world_t* world;
};

} /* namespace mcell */

#endif /* SRC4_DEFRAGMENTATION_EVENT_H_ */
