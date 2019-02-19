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

#ifndef SRC4_DIFFUSE_REACT_EVENT_H_
#define SRC4_DIFFUSE_REACT_EVENT_H_

#include "base_event.h"

namespace mcell {

// created in mcell3_world_converter::create_diffusion_events() before any other events,
// so even if release is created for time 0, this event is scheduled to occur as first one
class diffuse_react_event_t : public base_event_t {
	diffuse_react_event_t(world_t* world_, float_t diffusion_time_step_) :
		base_event_t(TYPE_DIFFUSE_REACT_EVENT),
		world(world_),
		diffusion_time_step(diffusion_time_step_) {

		// repeat this event each
		periodicity_interval = diffusion_time_step;
	}
	void step();
	void dump(const std::string indent) {
		// TODO
	}

	world_t* world;

	// this event diffuses all molecules that have this diffusion time_step
	float_t diffusion_time_step;
};

} // namespace mcell

#endif /* SRC4_DIFFUSE_REACT_EVENT_H_ */
