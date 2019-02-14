/*
 * event.h
 *
 *  Created on: Jan 30, 2019
 *      Author: adam
 */

#ifndef _EVENT_H_
#define _EVENT_H_


#include "defines.h"

namespace mcell {

/**
 * Base class for all events.
 */
class base_event_t {
public:
	base_event_t() :
		event_time(TIME_INVALID) {
	}
	virtual ~base_event_t() {};
	virtual void step() = 0;
	virtual void dump(const std::string ind) = 0;

	float_t event_time;
};

}

#endif /* V4SRC_EVENT_H_ */
