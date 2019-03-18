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

using namespace std;

namespace mcell {

world_t::world_t() {
	// TODO: initialize rest of members
	world_constants.partition_edge_length = PARTITION_EDGE_LENGTH_DEFAULT;
	world_constants.subpartitions_per_partition_dimension = SUBPARTITIONS_PER_PARTITION_DIMENSION_DEFAULT;
}

void world_t::init_simulation() {

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

	// just to make sure that we have an item for all the species
	for (species_t& s: species) {
		bimolecular_reactions_map.insert( std::make_pair(s.species_id, species_reaction_map_t()) );
	}
	assert(bimolecular_reactions_map.size() == species.size());

	world_constants.init(&bimolecular_reactions_map);

	cout << "Creating initial partition with " <<  world_constants.subpartitions_per_partition_dimension << "^3 subvolumes.";

	// create initial partition with center at 0,0,0 - we woud like to have the partitions all the same,
	// not depend on some random initialization
	uint32_t index = add_partition(vec3_t(0, 0, 0));
	assert(index == PARTITION_INDEX_INITIAL);


}


void world_t::dump() {
	world_constants.dump();
	// species
	species_t::dump_array(species);

}

static uint64_t determine_output_frequency(uint64_t iterations) {
	uint64_t frequency;

	if (iterations < 10)
		frequency = 1;
	else if (iterations < 1000)
		frequency = 10;
	else if (iterations < 100000)
		frequency = 100;
	else if (iterations < 10000000)
		frequency = 1000;
	else if (iterations < 1000000000)
		frequency = 10000;
	else
		frequency = 100000;

  return frequency;
}

bool world_t::run_simulation() {

	init_simulation(); // must be the first one

	dump();

	// create event that will terminate our simulation
	end_simulation_event_t* end_event = new end_simulation_event_t();
	end_event->event_time = iterations;
	scheduler.schedule_event(end_event);

	bool end_simulation = false;
	float_t time = TIME_SIMULATION_START;
	float_t previous_time;
	current_iteration = 0;
	uint32_t how_often_to_report = determine_output_frequency(iterations);

	cout << "Iterations: " << current_iteration << " of " << iterations << "\n";

	do {
		previous_time = time;
		time = scheduler.handle_next_event(end_simulation);

		if (time > previous_time) {
			current_iteration++;
			if (current_iteration % how_often_to_report == 0) {
				cout << "Iterations: " << current_iteration << " of " << iterations << "\n";
			}
		}
	} while (!end_simulation);

	cout << "Iteration " << current_iteration << ", simulation finished successfully\n";

	return true;
}

} /* namespace mcell */
