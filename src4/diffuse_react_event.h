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
#include <unordered_map>
#include <unordered_set>
#include "base_event.h"
#include "partition.h"
#include "reaction.h"

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
			partition_t* partition_ptr,
			const molecule_idx_t diffused_molecule_idx_,
			const molecule_idx_t colliding_molecule_idx_,
			reaction_t* rx_ptr,
			const float_t& time_,
			const vec3_t& position_)
		:
			partition(partition_ptr),
			diffused_molecule_idx(diffused_molecule_idx_),
			colliding_molecule_idx(colliding_molecule_idx_),
			rx(rx_ptr),
			time(time_),
			position(position_)
			{
	}

	partition_t* partition;
	molecule_idx_t diffused_molecule_idx;
	molecule_idx_t colliding_molecule_idx;
	reaction_t* rx;
	float_t time;
	vec3_t position;

  void dump(partition_t& p, const std::string ind) const;
  std::string to_string() const;
  static void dump_array(partition_t& p, const std::vector<molecules_collision_t>& vec);

};


class molecule_to_diffuse_t {
public:
	molecule_to_diffuse_t(
			const molecule_idx_t idx_,
			const float_t remaining_time_step_)
		:
			idx(idx_),
			remaining_time_step(remaining_time_step_)
			{
	}
	molecule_idx_t idx;
	float_t remaining_time_step;
};

// created in mcell3_world_converter::create_diffusion_events() before any other events,
// so even if release is created for time 0, this event is scheduled to occur as first one
class diffuse_react_event_t : public base_event_t {
public:
	diffuse_react_event_t(world_t* world_, const float_t diffusion_time_step_) :
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


	std::vector<molecule_to_diffuse_t> new_molecules_to_diffuse;

private:
	void diffuse_molecules(partition_t& p, std::vector< molecule_idx_t >& indices);
	void diffuse_single_molecule(partition_t& p, const molecule_idx_t vm, const float_t curr_time_step);
	void pick_displacement(float_t scale /*space step*/, vec3_t& displacement);
	void compute_displacement(species_t& sp, vec3_t& displacement, float_t remaining_time_step);

	ray_trace_state_t ray_trace(
			partition_t& p,
			volume_molecule_t& vm, // molecule that we are diffusing, we are changing its pos  and possibly also subvolume
			vec3_t& remaining_displacement, // in/out - recomputed if there was a reflection
			std::vector<molecules_collision_t>& molecule_collisions, // possible reactions in this part of way marching, ordered by time
			vec3_t& new_position,
			uint32_t& new_subpartition_index
	);

	bool collide_and_react_with_vol_mol(
			partition_t& p,
			molecules_collision_t& collision,
			vec3_t& displacement,
			float_t remaining_time_step
	);


	int test_bimolecular(
			reaction_t& rx,
			volume_molecule_t& a1,
			volume_molecule_t& a2
	);

	int outcome_bimolecular(
			partition_t& p,
			molecules_collision_t& collision,
			int path,
			float_t remaining_time_step
	);

	int outcome_products_random(
			partition_t& p,
			molecules_collision_t& collision,
			int path,
			float_t remaining_time_step
	);

	subpartition_mask_t& get_sp_species_reacting_mols_cached_data(
			uint32_t sp_index, volume_molecule_t& vm, partition_t& p);

	// cache for possible reactions, reset every time step() is called
	//std::unordered_map<uint32_t, std::unordered_map<species_id_t, subpartition_mask_t> > cache_sp_species_reacting_mols;
};

} // namespace mcell

#endif /* SRC4_DIFFUSE_REACT_EVENT_H_ */
