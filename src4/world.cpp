/*
 * world.cpp
 *
 *  Created on: Feb 10, 2019
 *      Author: adam
 */

#include "world.h"
#include "end_simulation_event.h"

namespace mcell {

bool world_t::run_simulation() {

	float_t time = TIME_SIMULATION_START;

	// create event that will terminate our simulation
	end_simulation_event_t* end_event = new end_simulation_event_t();
	end_event->event_time = time_unit * iterations; // TODO: check against MCELL - we need to execute the same nr of iterations, not one less
	scheduler.schedule_event(end_event);

	bool end_simulation = false;
	do {
		time = scheduler.handle_next_event(end_simulation);
	} while (!end_simulation);

	return true;
}

} /* namespace mcell */
