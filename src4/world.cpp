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
extern "C" {
#include "rng.h" // MCell 3
}

#include "world.h"
#include "end_simulation_event.h"

namespace mcell {

world_t::world_t() {
	// TODO: initialize rest of members
	world_constants.partition_edge_length = PARTITION_EDGE_LENGTH_DEFAULT;
	world_constants.subpartitions_per_partition_dimension = SUBPARTITIONS_PER_PARTITION_DIMENSION_DEFAULT;
	world_constants.init_subpartition_edge_length();
}

void world_t::init_simulation() {
	// create initial partition with center at 0,0,0 - we woud like to have the partitions all the same,
	// not depend on some random initialization
	uint32_t index = add_partition(vec3_t(0, 0, 0));
	assert(index == PARTITION_INDEX_INITIAL);

	// create map for fast reaction searches
	for (reaction_t& r: reactions) {
		assert(r.reactants.size() == 2); // only bimolecular reactions are supported now

		// TODO check - for now we only support one outcome of a reaction
		if (bimolecular_reactions_map.count(r.reactants[0].species_id) != 0) {
			assert(bimolecular_reactions_map[r.reactants[0].species_id].count(r.reactants[1].species_id) == 0);
		}

		bimolecular_reactions_map[r.reactants[0].species_id][r.reactants[1].species_id] = &r;
		bimolecular_reactions_map[r.reactants[1].species_id][r.reactants[0].species_id] = &r;
	}
}


void world_t::dump() {
	// species
	species_t::dump_array(species);
}

bool world_t::run_simulation() {

	init_simulation();

	dump();

	// create event that will terminate our simulation
	end_simulation_event_t* end_event = new end_simulation_event_t();
	end_event->event_time = world_constants.time_unit * iterations; // TODO: check against MCELL - we need to execute the same nr of iterations, not one less
	scheduler.schedule_event(end_event);

	bool end_simulation = false;
	float_t time = TIME_SIMULATION_START;
	float_t previous_time;
	uint32_t iteration = 0;
	do {
		previous_time = time;
		time = scheduler.handle_next_event(end_simulation);

		if (time > previous_time) {
			if (iteration % 100 == 0) {
				std::cout << "Iteration " << iteration << "\n";
			}
			iteration++;
		}
	} while (!end_simulation);

	std::cout << "Iteration " << iteration << ", simulation finished successfully\n";

	return true;
}

} /* namespace mcell */
