/*
 * end_simulation_event.h
 *
 *  Created on: Feb 14, 2019
 *      Author: adam
 */

#ifndef SRC4_END_SIMULATION_EVENT_H_
#define SRC4_END_SIMULATION_EVENT_H_

#include "base_event.h"

namespace mcell {

class end_simulation_event_t: public mcell::base_event_t {
public:
	end_simulation_event_t() :
		base_event_t(TYPE_END_SIMULATION_EVENT) {
	}
	void step() {
		// does nothing, this type of event is detected in the scheduler
	}
	void dump(const std::string indent) {
		// TODO
	}
};

} // namespace mcell

#endif /* SRC4_END_SIMULATION_EVENT_H_ */
