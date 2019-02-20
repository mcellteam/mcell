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

#ifndef SRC4_RELEASE_EVENT_H_
#define SRC4_RELEASE_EVENT_H_

#include "base_event.h"

namespace mcell {

class release_event_t: public base_event_t {
public:
	release_event_t(world_t* world_) :
		base_event_t(EVENT_TYPE_INDEX_RELEASE),
		species_id(SPECIES_ID_INVALID),
		release_number(0),
		world(world_) {
	}
	virtual ~release_event_t() {}

	virtual void step();
	virtual void dump(const std::string indent);

	vec3_t location;
	species_id_t species_id;
	uint32_t release_number; // number of molecules to release
	std::string name;

	world_t* world;
};

} // namespace mcell


#endif /* SRC4_RELEASE_EVENT_H_ */
