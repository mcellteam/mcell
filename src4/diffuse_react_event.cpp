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

#include <iostream>

#include "diffuse_react_event.h"
#include "world.h"
#include "partition.h"


using namespace std;

namespace mcell {

void diffuse_react_event_t::dump(const std::string indent) {
	cout << indent << "Diffuse-react event:\n";
	std::string ind2 = indent + "  ";
	base_event_t::dump(ind2);
	cout << ind2 << "diffusion_time_step: \t\t" << diffusion_time_step << " [float_t] \t\t\n";
}


void diffuse_react_event_t::step() {
	// for each partition
	for (partition_t& p: world->partitions) {
		// 1) select the right list of molecules from volume_molecule_indices_per_time_step
		uint32_t list_index = p.get_molecule_list_index_for_time_step(diffusion_time_step);
		if (list_index != PARTITION_INDEX_INVALID) {
			diffuse_molecules(p, p.volume_molecule_indices_per_time_step[list_index].second);
		}
	}
}


void diffuse_react_event_t::diffuse_molecules(partition_t& p, std::vector< molecule_index_t >& indices) {
	for (molecule_index_t i: indices) {
		volume_molecule_t& m = p.volume_molecules[i];
		if (m.is_defunct())
			continue;
		species_t& sp = world->species[m.species_id];

		// 2) diffuse each molecule
		// TBD: reflections
		vec3_t displacement;
		compute_displacement(sp, displacement);

		vec3_t new_pos = m.pos + displacement;

		// 3) detect collisions with other molecules
		// TODO

		// 4) evaluate and possible execute reactions
		// TODO

		// are we still in the same partition or do we need to move?
		bool move_to_another_partition = !p.in_this_partition(new_pos);
		assert(!move_to_another_partition && "TODO");

		// for now, and also to be able to compare with MCell 3,
		// let's move the particle right away
		// later, all particles will be moved after reactions were resolved
		m.pos = new_pos;
	}
}


void diffuse_react_event_t::pick_displacement(float_t scale /*space step*/, vec3_t& displacement) {
	displacement.x = scale * rng_gauss(&(world->rng)) * .70710678118654752440;
	displacement.y = scale * rng_gauss(&(world->rng)) * .70710678118654752440;
	displacement.z = scale * rng_gauss(&(world->rng)) * .70710678118654752440;
}


void diffuse_react_event_t::compute_displacement(species_t& sp, vec3_t& displacement) {
	pick_displacement(sp.space_step, displacement);
}



} /* namespace mcell */
