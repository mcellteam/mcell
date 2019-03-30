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

#ifndef SRC4_BASE_EVENT_H_
#define SRC4_BASE_EVENT_H_


#include "defines.h"

namespace mcell {

class world_t;

typedef int event_type_index_t;

// Value specifies ordering when two events are scheduled for exactly the same time
// The 'holes' are there on purpose for ordering of external events
const event_type_index_t EVENT_TYPE_INDEX_INVALID = -1;
const event_type_index_t EVENT_TYPE_INDEX_DIFFUSE_REACT = 100;
const event_type_index_t EVENT_TYPE_INDEX_RELEASE = 200;
const event_type_index_t EVENT_TYPE_INDEX_VIZ_OUTPUT = 300;
const event_type_index_t EVENT_TYPE_INDEX_DEFRAGMENTATION = 900;
const event_type_index_t EVENT_TYPE_INDEX_END_SIMULATION = 1000;

/**
 * Base class for all events.
 * Should be independent on the mcell world in order to
 * integrate also other simulators.
 */
class base_event_t {
public:
	base_event_t(const event_type_index_t t) :
		event_time(TIME_INVALID), periodicity_interval(0), type_index(t) {
	}
	virtual ~base_event_t() {};
	virtual void step() = 0;
	virtual void dump(const std::string indent);

	// time when this object;s step() method will be callled
	float_t event_time;

	// once this event is executed, schedule next one after this interval
	// do not schedule if the value is 0
	float_t periodicity_interval;

	// this value specifies both id of the event and ordering when multiple
	// events of a different type are sheduled for the same time
	event_type_index_t type_index;
};

}

#endif // SRC4_BASE_EVENT_H_
