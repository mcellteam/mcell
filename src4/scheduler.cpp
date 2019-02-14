/*
 * scheduler.cpp
 *
 *  Created on: Jan 30, 2019
 *      Author: adam
 */

#include "scheduler.h"

namespace mcell {

void calendar_t::insert(base_event_t* event) {



}

void scheduler_t::schedule_event(base_event_t* event) {
	calendar.insert(event);
}


} /* namespace mcell */
