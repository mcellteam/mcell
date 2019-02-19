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
