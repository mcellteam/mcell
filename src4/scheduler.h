/*
 * scheduler.h
 *
 *  Created on: Jan 30, 2019
 *      Author: adam
 */

#ifndef SRC4_SCHEDULER_H_
#define SRC4_SCHEDULER_H_

#include <deque>
#include <list>

#include "base_event.h"

namespace mcell {

const float_t BUCKET_TIME_INTERVAL = 1e-6; // time step - 1e-6

class calendar_t {
public:
	void insert(base_event_t* event);
	base_event_t* pop_next();

private:
	std::dequeue< std::list<base_event_t*>> > queue; // fix: bucket


};

class scheduler_t {
public:
	// scheduler becomes owner of the base_event object
	void schedule_event(base_event_t* event);


	calendar_t calendar;
};

} // namespace mcell

#endif /* SRC4_SCHEDULER_H_ */
