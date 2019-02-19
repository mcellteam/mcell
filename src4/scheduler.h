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

#ifndef SRC4_SCHEDULER_H_
#define SRC4_SCHEDULER_H_

#include <deque>
#include <list>

#include "base_event.h"

namespace mcell {

// TODO: what about rounding?
// prefferably, we should represent the time interval precisely
const float_t BUCKET_TIME_INTERVAL = 1e-6; // time step - 1e-6


class bucket_t {
public:
	bucket_t(float_t start_time_) :
		start_time(start_time_) {
	}
	~bucket_t();
	void insert(base_event_t* event);

	float_t start_time;
	std::list<base_event_t*> events;
};


typedef std::deque<bucket_t> bucket_deque_t;


class calendar_t {
public:
	calendar_t() {
		// create at least one item?
		queue.push_back(bucket_t(TIME_SIMULATION_START));
	}
	~calendar_t() {
		// implicitly calls destructors of items in queue and
		// deletes all events
	}

	void insert(base_event_t* event);
	base_event_t* pop_next();

private:
	float_t get_first_bucket_start_time() {
		assert(queue.size() != 0);
		return queue.front().start_time;
	}

	float_t event_time_to_bucket_start_time(const float_t time) {
		// flooring to a multiple of BUCKET_TIME_INTERVAL
		return (float_t)((int)(time / BUCKET_TIME_INTERVAL)) * BUCKET_TIME_INTERVAL;
	}

	bucket_deque_t::iterator get_or_create_bucket(const float_t time);

	// queue might be empty
	bucket_deque_t queue;
};


class scheduler_t {
public:
	// scheduler becomes owner of the base_event object
	void schedule_event(base_event_t* event);

	// returns new time
	float_t handle_next_event(bool &end_simulation);

	calendar_t calendar;
};

} // namespace mcell

#endif /* SRC4_SCHEDULER_H_ */
