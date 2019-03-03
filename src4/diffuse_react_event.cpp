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
		volume_molecule_t& vm = p.volume_molecules[i];
		if (vm.is_defunct())
			continue;
		species_t& spec = world->species[vm.species_id];


		// 2) diffuse each molecule - get information on position change
		// TBD: reflections
		vec3_t displacement;
		compute_displacement(spec, displacement);


		// ray_trace
		// - which subpartitions are crossed?


		// "housekeeping"

		// inertness?


	  /* scan subvolume for possible mol-mol reactions with vm */
		// for now, just collect all molecules
		// this is filtering out, do it later
/*		vector<collision_t> collisions;
	  if (spec.flags & SPECIES_FLAG_CAN_VOLVOL) {
	    determine_mol_mol_reactions(vm, p, collisions); <- this might belong to ray_trace
	  }*/


		// 3) detect collisions with other molecule
		//vec3_t orig_pos = vm.pos;
		vec3_t remaining_displacement = displacement;
		uint32_t new_subpartition_index;
		ray_trace_state_t state;
		std::vector<molecules_collision_t> molecule_collisions;
		do {
			state =
					ray_trace(
							p, vm /* changes position */,
							remaining_displacement,
							molecule_collisions, new_subpartition_index);

		} while (state != RAY_TRACE_FINISHED);

		// 4) evaluate and possible execute reactions
		// TODO

		// are we still in the same partition or do we need to move?
		bool move_to_another_partition = !p.in_this_partition(vm.pos);
		assert(!move_to_another_partition && "TODO");

		// change subpartition
		p.change_molecule_subpartition(vm, new_subpartition_index);

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

/*
bool ray_trace_subpartition_boundary(
		volume_molecule& vm, // molecule that we are diffusing, we are changing its pos  and possibly also subvolume
		vec3_t& remaining_displacement, // in/out - recomputed if there was a reflection

)
*/

// collect subpartition indices that we are crossing
void diffuse_react_event_t::collect_crossed_subpartitions(
	const partition_t& p,
	volume_molecule_t& vm, // molecule that we are diffusing, we are changing its pos  and possibly also subvolume
	vec3_t& displacement, // in/out - recomputed if there was a reflection
	std::vector<uint32_t>& crossed_subparition_indices) {

	crossed_subparition_indices.clear();
	// remeber the starting subpartition
	crossed_subparition_indices.push_back(vm.subpartition_index);

	// destination
	vec3_t dest_pos = vm.pos + displacement;

	// get 3d indices of start and end subpartitions
	ivec3_t src_sp_indices, dest_sp_indices;
	p.get_subpartition_3d_indices_from_index(vm.subpartition_index, src_sp_indices);
	p.get_subpartition_3d_indices(dest_pos, dest_sp_indices);

	// collect subpartitions on the way
	ivec3_t sp_indices_abs = glm::abs(dest_sp_indices - src_sp_indices);
	int sp_indices_sum = sp_indices_abs.x + sp_indices_abs.y + sp_indices_abs.z;
	if (sp_indices_sum > 0) {

		uint32_t dest_sp_index = p.get_subpartition_index_from_3d_indices(dest_sp_indices);
		uint32_t curr_sp_index;

		if (sp_indices_sum == 1) {
			// crossing just one boundary, simply add destination index
			crossed_subparition_indices.push_back(dest_sp_index);
		}
		else {
			// directions - 0/1 for multiplication
			ivec3_t dir_urb_multiplier = ivec3_t(glm::greaterThan(displacement, vec3_t(0)));
			assert(dir_urb_multiplier.x == 0 || dir_urb_multiplier.x == 1);
			assert(dir_urb_multiplier.y == 0 || dir_urb_multiplier.y == 1);
			assert(dir_urb_multiplier.z == 0 || dir_urb_multiplier.z == 1);

			ivec3_t dir_urb_addend;
			dir_urb_addend.x = (dir_urb_multiplier.x == 0) ? -1 : 1;
			dir_urb_addend.y = (dir_urb_multiplier.y == 0) ? -1 : 1;
			dir_urb_addend.z = (dir_urb_multiplier.z == 0) ? -1 : 1;

			vec3_t curr_pos = vm.pos;
			ivec3_t curr_sp_indices = src_sp_indices;

			do {
				// subpartition edges
				// = origin + subparition index * length + is_urb * length
				vec3_t sp_edges =
						p.origin_corner
						+	vec3_t(curr_sp_indices) * vec3_t(world->world_constants.subpartition_edge_length) // llf edge
						+ vec3_t(dir_urb_multiplier) * vec3_t(world->world_constants.subpartition_edge_length); // move if we go urb

				// compute time for the next subpartition collision, let's assume that displacemnt
				// is our speed vector and the total time to travel is 1
				//
				// pos(time) = pos + displacement * time, therefore
				// time = (pos(time) - vm.pos) / displacement
				// =>
				// time_to_subpart_edge = (subpart_edge - vm.pos) / displacement_speed
				vec3_t coll_times =
						(sp_edges - curr_pos) / displacement; // TODO: what if displacement is 0

				// which of the times is the smallest? - i.e. which boundary we hit first
				if (coll_times.x < coll_times.y && coll_times.x <= coll_times.z) {
					// new position on the edge of the subpartition
					curr_pos += displacement * coll_times.x;
					// and also update the xyz subpartition index
					curr_sp_indices.x += dir_urb_addend.x;
				}
				else if (coll_times.y <= coll_times.z) {
					// y
					curr_pos += displacement * coll_times.y;
					curr_sp_indices.y += dir_urb_addend.y;
				}
				else {
					// z
					curr_pos += displacement * coll_times.z;
					curr_sp_indices.z += dir_urb_addend.z;
				}

				curr_sp_index = p.get_subpartition_index_from_3d_indices(curr_sp_indices);

				crossed_subparition_indices.push_back(curr_sp_index);
			} while (curr_sp_index != dest_sp_index);
		}
	}
}


/***************************************************************************
collide_mol:
  In: starting coordinate
      vector to move along
      molecule we're checking for a collision
      double to store time of collision
      vector to store the location of the collision
  Out: True if collision was detected
  Note: collision_time and/or collision_pos may be modified even if there is no collision
        Not highly optimized yet.
***************************************************************************/
bool diffuse_react_event_t::collide_mol(
		volume_molecule_t& diffused_vm,
		vec3_t& move,
    volume_molecule_t& colliding_vm, //a
		float_t& rel_collision_time, //t
		vec3_t& rel_collision_pos, // ???
    float_t rx_radius_3d) {

  //if ((a->properties->flags & ON_GRID) != 0)
  //  return COLLIDE_MISS; /* Should never call on surface molecule! */

  vec3_t& pos = colliding_vm.pos; /* Position of target molecule */
  vec3_t dir = pos - diffused_vm.pos;  /* From starting point of moving molecule to target */

  float_t d = glm::dot((glm_vec3_t)dir, (glm_vec3_t)move);        /* Dot product of movement vector and vector to target */
  float_t sigma2 = rx_radius_3d * rx_radius_3d;   /* Square of interaction radius */

  /* Miss the molecule if it's behind us */
  if (d < 0)
    return false;

  float_t movelen2 = glm::dot((glm_vec3_t)move, (glm_vec3_t)move); /* Square of distance the moving molecule travels */

  /* check whether the test molecule is further than the displacement. */
  if (d > movelen2)
    return false;

  float_t dirlen2 = glm::dot((glm_vec3_t)dir, (glm_vec3_t)dir);

  /* check whether the moving molecule will miss interaction disk of the
     test molecule.*/
  if (movelen2 * dirlen2 - d * d > movelen2 * sigma2)
    return false;

  rel_collision_time = d / movelen2; // not exact time, only for comparision, exact time should be d/sqrt(movelen2*dirlen2);

  rel_collision_pos = diffused_vm.pos + rel_collision_time * move;
  return COLLIDE_VOL_M;
}

// collect possible collisions until a wall is hit
ray_trace_state_t diffuse_react_event_t::ray_trace(
		partition_t& p,
		volume_molecule_t& vm, // molecule that we are diffusing, we are changing its pos  and possibly also subvolume
		vec3_t& remaining_displacement, // in/out - recomputed if there was a reflection
		std::vector<molecules_collision_t>& molecule_collisions, // possible reactions in this part of way marching, ordered by time
		uint32_t& new_subpartition_index
		) {

  // first get what subpartitions might be relevant
	std::vector<uint32_t> crossed_subparition_indices;
	collect_crossed_subpartitions(p, vm, remaining_displacement, crossed_subparition_indices);


	// TBD: check wall collisions
	// here we can return RAY_TRACE_HIT_WALL

	// for each SP
	for (uint32_t sp_index: crossed_subparition_indices) {
		subpartition_mask_t& sp_mask = p.volume_molecules_subpartition_masks[sp_index];
		// for each molecule in this SP
		for (uint32_t vm_index: sp_mask) {

			volume_molecule_t& colliding_vm = p.volume_molecules[vm_index];
			if (colliding_vm.is_defunct()) {
				continue;
			}

			// can we react?
			// TODO: filter out molecules that cannot react right away

			// we would like to compute everything that's needed just once
			float_t time;
			vec3_t position;
			// collide_mol must be inlined because many things are computed all over there
			if (collide_mol(vm, remaining_displacement, colliding_vm, time, position, world->world_constants.rx_radius_3d)) {
				molecule_collisions.push_back(
						molecules_collision_t(vm, colliding_vm, time, position)
				);
			}
		}
	}

	// the value is valid only when RAY_TRACE_FINISHED is returned
	new_subpartition_index = crossed_subparition_indices.back();
	vm.pos = vm.pos + remaining_displacement;

  return RAY_TRACE_FINISHED; // no wall was hit
}


#if 0
int diffuse_react_event_t::trigger_bimolecular(
		volume_molecule_t& diffused_mol,
		volume_molecule_t& colliding_mol,
		std::vector<molecules_collision_t>& possible_collisions
) {
  //int num_matching_rxns = 0; /* number of matching reactions */

  // seems to only check whether there can be a reaction
  // TODO: for now, we consider all molecules
  possible_collisions.push_back(molecules_collision_t(colliding_mol));
  return 1;
}

/******************************************************************************
 * the determine_mol_mol_reactions helper function is used in diffuse_3D to
 * compute all possible molecule molecule reactions between the diffusing
 * molecule m and all other volume molecules in the subvolume.
 ******************************************************************************/
void diffuse_react_event_t::determine_mol_mol_reactions(
		volume_molecule_t& diffused_mol,
		partition_t& p,
		std::vector<molecules_collision_t>& possible_collisions) {

	// note: this for loop goes over species in MCell 3
	// go through all items (volume molecules for now) in a subpartition
	subpartition_mask_t& volume_molecules_mask =
			p.volume_molecules_subpartition_masks[diffused_mol.subpartition_index];

	for (uint32_t i = 0; i < volume_molecules_mask.size(); i++) {
		// FIXME: this is a not very efficient way how to iterate through a bitset
		// need a different bitset implementation
		if (!volume_molecules_mask.test(i)) {
			continue;
		}
		volume_molecule_t& colliding_mol = p.volume_molecules[i];

		if (colliding_mol.is_defunct()) {
			continue;
		}

		// same molecule?
		if (&diffused_mol == &colliding_mol) {
			continue;
		}

		// are there possible reactions between these species?
		// TODO: cache this result somewhere
		if (!world->can_react(diffused_mol.species_id, colliding_mol.species_id)) {
			continue;
		}

		// check collisions
		trigger_bimolecular(diffused_mol, colliding_mol, possible_collisions);
	}

}
#endif

} /* namespace mcell */
