/*
 * scheduler.cpp
 *
 *  Created on: Jan 30, 2019
 *      Author: adam
 */

#include "scheduler.h"

namespace mcell {



void bucket_t::insert(base_event_t* event) {
	// check right away if it does not belong to the end
	if (events.empty() || events.back()->event_time <= event->event_time) {
		events.push_back(event);
	}
	else {
		// simply go through the list and put our item to the right time
		auto it = events.begin();
		while (it != events.end() && (*it)->event_time <= event->event_time) {
			it++;
		}
		events.insert(it, event);
	}
}

bucket_t::~bucket_t() {
  for (auto it = events.begin(); it != events.end(); it++) {
    delete *it;
  }
}


void calendar_t::insert(base_event_t* event) {
	float_t bucket_start_time = event_time_to_bucket_start_time(event->event_time);
	if (queue.empty()) {
		// no items yet - simply create new bucket and insert our event there
		queue.push_back( bucket_t(bucket_start_time) );
		queue.front().insert(event);
	}
	else {
		// we first need to find out whether we already have a bucket for this event
		float_t first_start_time = get_first_bucket_start_time();
		assert(bucket_start_time - first_start_time >= 0 && "cannot schedule to the past"); // some eps?
		size_t buckets_from_first = (bucket_start_time - first_start_time) / BUCKET_TIME_INTERVAL;

		if (buckets_from_first < queue.size()) {
			// bucket exists
			queue[buckets_from_first].insert(event);
		}
		else {
			// we need to create new buckets
			// note: MCell 3 uses logarithmic scale
			size_t missing_buckets = buckets_from_first - queue.size() + 1;
			float_t next_time = queue.back().start_time + BUCKET_TIME_INTERVAL;
			for (size_t i = 0; i < missing_buckets; i++) {
				queue.push_back(bucket_t(next_time));
				next_time += BUCKET_TIME_INTERVAL;
			}
		}

	}
}


base_event_t* calendar_t::pop_next() {
	while (queue.front().events.empty()) {
		queue.pop_front();
	}
	base_event_t* next_event = queue.front().events.front();
	queue.front().events.pop_front();
	return next_event;
}


void scheduler_t::schedule_event(base_event_t* event) {
	calendar.insert(event);
}


float_t scheduler_t::handle_next_event(bool &end_simulation) {

	base_event_t* event = calendar.pop_next();
	assert(event != NULL && "Empty event queue - at least end simulation event should be present");
	float_t event_time = event->event_time;

	event->step();

	// schedule itself for the next period or just delete
	if (event->periodicity_interval != 0) {
		event->event_time += event->periodicity_interval;
		calendar.insert(event);
	}
	else {
		delete event;
	}

	end_simulation = event->type == TYPE_END_SIMULATION_EVENT;
	return event_time;
}


} /* namespace mcell */
