/******************************************************************************
 *
 * Copyright (C) 2019 by
 * The Salk Institute for Biological Studies and
 * Pittsburgh Supercomputing Center, Carnegie Mellon University
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
 * USA.
 *
******************************************************************************/

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
