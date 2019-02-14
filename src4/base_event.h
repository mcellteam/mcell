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

class world_t;

enum event_type_t {
	TYPE_INVALID_EVENT,
	TYPE_END_SIMULATION_EVENT,
	TYPE_RELEASE_EVENT,
	TYPE_DIFFUSE_REACT_EVENT
};

/**
 * Base class for all events.
 * Should be independent on the mcell world in order to
 * integrate also other simulators.
 */
class base_event_t {
public:
	base_event_t(const event_type_t t) :
		event_time(TIME_INVALID), periodicity_interval(0), type(t) {
	}
	virtual ~base_event_t() {};
	virtual void step() = 0;
	virtual void dump(const std::string indent) = 0;

	float_t event_time;

	// once this event is executed, schedule next one after this interval
	// do not schedule if the value is 0
	float_t periodicity_interval;

	// we do not want to use C++ RTTI (runtime type identification),
	// instead each superclass sets its own type
	event_type_t type;

};

}

#endif /* V4SRC_EVENT_H_ */
