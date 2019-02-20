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

#include <iostream>

#include "release_event.h"
#include "world.h"
#include "partition.h"

using namespace std;

namespace mcell {

void release_event_t::dump(const std::string ind) {
	cout << "Release event:\n";
	std::string ind2 = ind + "  ";
	base_event_t::dump(ind2);
	cout << ind2 << "location: \t\t" << location << " [vec3_t] \t\t\n";
	cout << ind2 << "species_id: \t\t" << species_id << " [species_id_t] \t\t\n";
	cout << ind2 << "release_number: \t\t" << release_number << " [uint32_t] \t\t\n";
	cout << ind2 << "name: \t\t" << name << " [string] \t\t\n";
}


void release_event_t::step() {
	// for now, let's simply release 'release_number' of molecules of 'species_id'
	// at 'location'

	partition_t& p = world->partitions[world->get_or_add_partition_index(location)];
	volume_molecule_t m(species_id, location);
	float_t time_step = world->species[species_id].time_step;
	uint32_t time_step_index = p.get_or_add_molecule_list_index_for_time_step(time_step);

	for (uint32_t i = 0; i < release_number; i++) {
		p.add_volume_molecule(m, time_step_index);
	}

}

} // namespace mcell


