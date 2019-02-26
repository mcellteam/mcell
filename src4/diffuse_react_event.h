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

#include <vector>
#include "base_event.h"


namespace mcell {

class partition_t;
class volume_molecule_t;
class species_t;


class collision_t {
public:
	collision_t(volume_molecule_t& colliding_molecule_)
		: colliding_molecule(colliding_molecule_) {
	}
	volume_molecule_t& colliding_molecule;
};

// created in mcell3_world_converter::create_diffusion_events() before any other events,
// so even if release is created for time 0, this event is scheduled to occur as first one
class diffuse_react_event_t : public base_event_t {
public:
	diffuse_react_event_t(world_t* world_, float_t diffusion_time_step_) :
		base_event_t(EVENT_TYPE_INDEX_DIFFUSE_REACT),
		world(world_),
		diffusion_time_step(diffusion_time_step_) {

		// repeat this event each
		periodicity_interval = diffusion_time_step;
	}
	void step();
	void dump(const std::string indent);

	world_t* world;

	// this event diffuses all molecules that have this diffusion time_step
	float_t diffusion_time_step;

private:
	void diffuse_molecules(partition_t& p, std::vector< molecule_index_t >& indices);
	void pick_displacement(float_t scale /*space step*/, vec3_t& displacement);
	void compute_displacement(species_t& sp, vec3_t& displacement);

	int trigger_bimolecular(
			volume_molecule_t& diffused_mol,
			volume_molecule_t& colliding_mol,
			std::vector<collision_t>& possible_collisions);

	void determine_mol_mol_reactions(
			volume_molecule_t& vm,
			partition_t& p,
			std::vector<collision_t>& possible_collisions);

	void ray_trace(
			volume_molecule& diffused_mol,
			vec3_t& displacement,
			std::vector<collision_t>& possible_collisions);

};

} // namespace mcell

#endif /* SRC4_DIFFUSE_REACT_EVENT_H_ */
