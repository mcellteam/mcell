/*
 * diffuse.h
 *
 *  Created on: Jan 30, 2019
 *      Author: adam
 */

#ifndef SRC4_DIFFUSE_REACT_EVENT_H_
#define SRC4_DIFFUSE_REACT_EVENT_H_

#include "base_event.h"

namespace mcell {

class diffuse_react_event_t : public base_event_t {
	diffuse_react_event_t(world_t* world_) :
		base_event_t(TYPE_DIFFUSE_REACT_EVENT),
		world(world_) {
	}
	void step() {
		// TODO
	}
	void dump(const std::string indent) {
		// TODO
	}

	world_t* world;
};

} // namespace mcell

#endif /* SRC4_DIFFUSE_REACT_EVENT_H_ */
