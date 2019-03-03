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


enum ray_trace_state_t {
	RAY_TRACE_HIT_UNDEFINED,
	RAY_TRACE_HIT_SUBPARTITION,
	RAY_TRACE_HIT_WALL,
	RAY_TRACE_FINISHED
};


class molecules_collision_t {
public:
	molecules_collision_t(
			volume_molecule_t& diffused_molecule_ref,
			volume_molecule_t& colliding_molecule_ref,
			float_t& time_,
			vec3_t& position_)
		: diffused_molecule(diffused_molecule_ref),
			colliding_molecule(colliding_molecule_ref),
			time(time_),
			position(position_)
			{
	}
	volume_molecule_t& diffused_molecule;
	volume_molecule_t& colliding_molecule;
	float_t time;
	vec3_t position;
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

	/*int trigger_bimolecular(
			volume_molecule_t& diffused_mol,
			volume_molecule_t& colliding_mol,
			std::vector<molecules_collision_t>& possible_collisions);

	void determine_mol_mol_reactions(
			volume_molecule_t& vm,
			partition_t& p,
			std::vector<molecules_collision_t>& possible_collisions);
*/
	void collect_crossed_subpartitions(
		const partition_t& p,
		volume_molecule_t& vm, // molecule that we are diffusing, we are changing its pos  and possibly also subvolume
		vec3_t& remaining_displacement, // in/out - recomputed if there was a reflection
		std::vector<uint32_t>& crossed_subparition_indices
	);

	bool collide_mol(
			volume_molecule_t& diffused_vm,
			vec3_t& move,
	    volume_molecule_t& colliding_vm,
			float_t& rel_collision_time,
			vec3_t& rel_collision_pos,
	    float_t rx_radius_3d
	);

	ray_trace_state_t ray_trace(
			partition_t& p,
			volume_molecule_t& vm, // molecule that we are diffusing, we are changing its pos  and possibly also subvolume
			vec3_t& remaining_displacement, // in/out - recomputed if there was a reflection
			std::vector<molecules_collision_t>& molecule_collisions, // possible reactions in this part of way marching, ordered by time
			uint32_t& new_subpartition_index
	);

};

} // namespace mcell

#endif /* SRC4_DIFFUSE_REACT_EVENT_H_ */
