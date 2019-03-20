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
#include "mcell_structs.h"
}

#include <iostream>
#include <sstream>
#include <algorithm>

#include "diffuse_react_event.h"
#include "world.h"
#include "partition.h"

#include "debug_config.h"

using namespace std;

namespace mcell {

void diffuse_react_event_t::dump(const std::string indent) {
	cout << indent << "Diffuse-react event:\n";
	std::string ind2 = indent + "  ";
	base_event_t::dump(ind2);
	cout << ind2 << "diffusion_time_step: \t\t" << diffusion_time_step << " [float_t] \t\t\n";
}


void diffuse_react_event_t::step() {
	assert(world->partitions.size() == 1 && "Must extend cache to handle multiple partitions");

	// for each partition
	for (partition_t& p: world->partitions) {
		// 1) select the right list of molecules from volume_molecule_indices_per_time_step
		uint32_t list_index = p.get_molecule_list_index_for_time_step(diffusion_time_step);
		if (list_index != PARTITION_INDEX_INVALID) {
			diffuse_molecules(p, p.volume_molecule_indices_per_time_step[list_index].second);
		}
	}
}


void diffuse_react_event_t::diffuse_molecules(partition_t& p, std::vector< molecule_idx_t >& indices) {

	// Possible optimization:
	// GPU version - can make sets of possible collisions for a given species before
	// the whole diffuse is run
	// ! CPU version - compute sets of possible colilisions for a given species and "cache" it

	// first diffuse already existing molecules
	uint32_t existing_mols_count = indices.size();
	for (uint32_t i = 0; i < existing_mols_count; i++) {
		molecule_idx_t idx = indices[i];
		// existing molecules - simulate whole time step
		diffuse_single_molecule(p, idx, diffusion_time_step);
	}

	// need to call .size() each iteration because the size can increase
	for (uint32_t i = 0; i < new_molecules_to_diffuse.size(); i++) {
		// new molecules created with reactions - simulate remaining time
		diffuse_single_molecule(p, new_molecules_to_diffuse[i].idx, new_molecules_to_diffuse[i].remaining_time_step);
	}

	new_molecules_to_diffuse.clear();
}

void diffuse_react_event_t::diffuse_single_molecule(partition_t& p, const molecule_idx_t vm_idx, const float_t remaining_time_step) {

	volume_molecule_t& vm = p.volume_molecules[vm_idx];

	if (vm.is_defunct())
		return;

	species_t& spec = world->species[vm.species_id];

#ifdef DEBUG_DIFFUSION
	DUMP_CONDITION4(
		vm.dump(world, "", "Diffusing vm:", world->current_iteration);
	);
#endif

	// diffuse each molecule - get information on position change
	// TBD: reflections
	vec3_t displacement;
	compute_displacement(spec, displacement, remaining_time_step);

#ifdef DEBUG_DIFFUSION
	DUMP_CONDITION4(
		displacement.dump("  displacement:", "");
	);
#endif

	// 3) detect collisions with other molecule
	//vec3_t orig_pos = vm.pos;
	vec3_t remaining_displacement = displacement;
	uint32_t new_subpartition_index;
	vec3_t new_position;
	ray_trace_state_t state;
	std::vector<molecules_collision_t> molecule_collisions;
	bool was_defunct = false;
	do {
		state =
				ray_trace(
						p, vm /* changes position */,
						remaining_displacement,
						molecule_collisions,
						new_position,
						new_subpartition_index
				);

		// sort collisions by time
		// note: it will suffice to sort them once, this is here mainly because of the dump
		std::sort( molecule_collisions.begin(), molecule_collisions.end(), [ ]( const auto& lhs, const auto& rhs )
		{
		  return lhs.time < rhs.time;
		});

#ifdef DEBUG_COLLISIONS
		DUMP_CONDITION4(
			molecules_collision_t::dump_array(p, molecule_collisions);
		);
#endif

	} while (state != RAY_TRACE_FINISHED);

	// 4) evaluate and possible execute reactions
	// they are supposed to be sorted in reverse by time
	//for (int collision_idx = molecule_collisions.size() - 1; collision_idx >= 0; collision_idx--) {
	for (size_t collision_idx = 0; collision_idx < molecule_collisions.size(); collision_idx++) {
		molecules_collision_t& collision = molecule_collisions[collision_idx];

		assert(collision.time > 0 && collision.time <= 1);

		// mol-mol collision
		if (collision.time < EPS) {
				continue;
		}

		// for now. do the change right away, but we will need to
		// cache these changes and do them right away
		// might invalidate references!
		if (collide_and_react_with_vol_mol(p, collision, displacement, remaining_time_step)) {
			// molecule was destroyed
   		was_defunct = true;
			break;
		}

	}

	if (!was_defunct) {
		// need to get a new reference
		volume_molecule_t& vm_refresh = p.volume_molecules[vm_idx];

		// finally move molecule to its destination
		vm_refresh.pos = new_position;

		// are we still in the same partition or do we need to move?
		bool move_to_another_partition = !p.in_this_partition(vm_refresh.pos);
		assert(!move_to_another_partition && "TODO");

		// change subpartition
		p.change_molecule_subpartition(vm_refresh, new_subpartition_index);
	}
}

void diffuse_react_event_t::pick_displacement(float_t scale /*space step*/, vec3_t& displacement) {
	displacement.x = scale * rng_gauss(&(world->rng)) * .70710678118654752440;
	displacement.y = scale * rng_gauss(&(world->rng)) * .70710678118654752440;
	displacement.z = scale * rng_gauss(&(world->rng)) * .70710678118654752440;
}


void diffuse_react_event_t::compute_displacement(species_t& sp, vec3_t& displacement, float_t remaining_time_step) {
  float_t rate_factor = (remaining_time_step == 1.0) ? 1.0 : sqrt(remaining_time_step);
	pick_displacement(sp.space_step * rate_factor, displacement);
}


void diffuse_react_event_t::collect_neigboring_subparitions(
		const partition_t& p,
		vec3_t& pos,
		ivec3_t& sp_indices,
		std::set<uint32_t>& crossed_subparition_indices
) {
	float_t r = world->world_constants.rx_radius_3d;
	float_t sp_len = world->world_constants.subpartition_edge_length;

	vec3_t rel_pos = pos - p.origin_corner;

	// check neighbors
	// TODO: boundaries must be precomputed - would ireally help? it's just a few multiplications
	// x subpart boundaries
	// left (x)
	int x_dir_used = 0;
	float_t x_boundary = sp_indices.x * sp_len;
	if (rel_pos.x - r < x_boundary) {
		crossed_subparition_indices.insert(
				p.get_subpartition_index_from_3d_indices(sp_indices.x - 1, sp_indices.y, sp_indices.z));
		x_dir_used = -1;
	}
	// right (x)
	else if (rel_pos.x + r > x_boundary + sp_len) { // assuming that subpartitions are larger than radius
		crossed_subparition_indices.insert(
				p.get_subpartition_index_from_3d_indices(sp_indices.x + 1, sp_indices.y, sp_indices.z));
		x_dir_used = +1;
	}

	// upper (y)
	int y_dir_used = 0;
	float_t y_boundary = sp_indices.y * sp_len;
	if (rel_pos.y - r < y_boundary) {
		crossed_subparition_indices.insert(
				p.get_subpartition_index_from_3d_indices(sp_indices.x, sp_indices.y - 1, sp_indices.z));
		y_dir_used = -1;
	}
	// right (y)
	else if (rel_pos.y + r > y_boundary + sp_len) {
		crossed_subparition_indices.insert(
				p.get_subpartition_index_from_3d_indices(sp_indices.x, sp_indices.y + 1, sp_indices.z));
		y_dir_used = +1;
	}

	// front (z)
	int z_dir_used = 0;
	float_t z_boundary = sp_indices.z * sp_len;
	if (rel_pos.z - r < z_boundary) {
		crossed_subparition_indices.insert(
				p.get_subpartition_index_from_3d_indices(sp_indices.x, sp_indices.y, sp_indices.z - 1));
		z_dir_used = -1;
	}
	// back (z)
	else if (rel_pos.z + r > z_boundary + sp_len) {
		crossed_subparition_indices.insert(
				p.get_subpartition_index_from_3d_indices(sp_indices.x, sp_indices.y, sp_indices.z + 1));
		z_dir_used = +1;
	}

	// we also have to count with movement in multiple dimensions

	// xy
	if (x_dir_used != 0 && y_dir_used != 0) {
		crossed_subparition_indices.insert(
				p.get_subpartition_index_from_3d_indices(sp_indices.x + x_dir_used, sp_indices.y + y_dir_used, sp_indices.z));
	}

	// xz
	if (x_dir_used != 0 && z_dir_used != 0) {
		crossed_subparition_indices.insert(
				p.get_subpartition_index_from_3d_indices(sp_indices.x + x_dir_used, sp_indices.y, sp_indices.z + z_dir_used));
	}

	// yz
	if (y_dir_used != 0 && z_dir_used != 0) {
		crossed_subparition_indices.insert(
				p.get_subpartition_index_from_3d_indices(sp_indices.x, sp_indices.y + y_dir_used, sp_indices.z + z_dir_used));
	}

	// xyz
	if (x_dir_used != 0 && y_dir_used != 0 && z_dir_used != 0) {
		crossed_subparition_indices.insert(
				p.get_subpartition_index_from_3d_indices(sp_indices.x + x_dir_used, sp_indices.y + y_dir_used, sp_indices.z + z_dir_used));
	}
}


// collect subpartition indices that we are crossing
void diffuse_react_event_t::collect_crossed_subpartitions(
	const partition_t& p,
	volume_molecule_t& vm, // molecule that we are diffusing, we are changing its pos  and possibly also subvolume
	vec3_t& displacement, // in/out - recomputed if there was a reflection
	std::set<uint32_t>& crossed_subparition_indices,
	uint32_t& last_subpartition_index
) {

	crossed_subparition_indices.clear();
	// remeber the starting subpartition
	crossed_subparition_indices.insert(vm.subpartition_index);

	// destination
	vec3_t dest_pos = vm.pos + displacement;

	// urb - upper, right, bottom
	ivec3_t dir_urb_direction = ivec3_t(glm::greaterThan(displacement, vec3_t(0)));
	assert(dir_urb_direction.x == 0 || dir_urb_direction.x == 1);
	assert(dir_urb_direction.y == 0 || dir_urb_direction.y == 1);
	assert(dir_urb_direction.z == 0 || dir_urb_direction.z == 1);

	// get 3d indices of start and end subpartitions
	ivec3_t src_sp_indices, dest_sp_indices;
	p.get_subpartition_3d_indices_from_index(vm.subpartition_index, src_sp_indices);
	p.get_subpartition_3d_indices(dest_pos, dest_sp_indices);

	// first check what's around the starting point
	collect_neigboring_subparitions(p, vm.pos, src_sp_indices, crossed_subparition_indices);

	// collect subpartitions on the way by always finding the point where a subpartition boundary is hit
	// we must do it eve when we are crossing just one subpartition because we might hit others while
	// moving along them
	if ( !glm::all( glm::equal(dest_sp_indices, src_sp_indices) ) ) {

		uint32_t dest_sp_index = p.get_subpartition_index_from_3d_indices(dest_sp_indices);
		last_subpartition_index = dest_sp_index;

		ivec3_t dir_urb_addend;
		dir_urb_addend.x = (dir_urb_direction.x == 0) ? -1 : 1;
		dir_urb_addend.y = (dir_urb_direction.y == 0) ? -1 : 1;
		dir_urb_addend.z = (dir_urb_direction.z == 0) ? -1 : 1;

		vec3_t curr_pos = vm.pos;
		ivec3_t curr_sp_indices = src_sp_indices;

		uint32_t curr_sp_index;

		vec3_t displacement_rcp = 1.0/displacement; // TODO: what if displacement is 0

		do {
			// subpartition edges
			// = origin + subparition index * length + is_urb * length
			vec3_t sp_edges =
					p.origin_corner
					+	vec3_t(curr_sp_indices) * vec3_t(world->world_constants.subpartition_edge_length) // llf edge
					+ vec3_t(dir_urb_direction) * vec3_t(world->world_constants.subpartition_edge_length); // move if we go urb

			// compute time for the next subpartition collision, let's assume that displacemnt
			// is our speed vector and the total time to travel is 1
			//
			// pos(time) = pos + displacement * time, therefore
			// time = (pos(time) - vm.pos) / displacement
			// =>
			// time_to_subpart_edge = (subpart_edge - vm.pos) / displacement_speed
			assert(displacement.x != 0 && displacement.y != 0 && displacement.z != 0);
			vec3_t coll_times = (sp_edges - curr_pos) * displacement_rcp;

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

			crossed_subparition_indices.insert(curr_sp_index);

			// also neighbors
			collect_neigboring_subparitions(p, curr_pos, curr_sp_indices, crossed_subparition_indices);

		} while (curr_sp_index != dest_sp_index);
	}
	else {
		// subpartition index did not change
		last_subpartition_index = vm.subpartition_index;
	}

	// finally check also neighbors in destination
	collect_neigboring_subparitions(p, dest_pos, dest_sp_indices, crossed_subparition_indices);
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

  /* Miss the molecule if it's behind us */
  if (d < 0) {
    return false;
  }

  float_t movelen2 = glm::dot((glm_vec3_t)move, (glm_vec3_t)move); /* Square of distance the moving molecule travels */

  /* check whether the test molecule is further than the displacement. */
  if (d > movelen2) {
    return false;
  }

  /* check whether the moving molecule will miss interaction disk of the
     test molecule.*/
  float_t dirlen2 = glm::dot((glm_vec3_t)dir, (glm_vec3_t)dir);
  float_t sigma2 = rx_radius_3d * rx_radius_3d;   /* Square of interaction radius */
  if (movelen2 * dirlen2 - d * d > movelen2 * sigma2) {
    return false;
  }

  /* reject collisions with itself */
  if (diffused_vm.idx == colliding_vm.idx) {
  	return false;
  }

  /* defunct - not probable */
	if (colliding_vm.is_defunct()) {
  	return false;
	}

  rel_collision_time = d / movelen2;

  rel_collision_pos = diffused_vm.pos + rel_collision_time * move;
  return COLLIDE_VOL_M;
}

// collect possible collisions until a wall is hit
ray_trace_state_t diffuse_react_event_t::ray_trace(
		partition_t& p,
		volume_molecule_t& vm, // molecule that we are diffusing, we are changing its pos  and possibly also subvolume
		vec3_t& remaining_displacement, // in/out - recomputed if there was a reflection
		std::vector<molecules_collision_t>& molecule_collisions, // possible reactions in this part of way marching, ordered by time
		vec3_t& new_position,
		uint32_t& new_subpartition_index
		) {

  // first get what subpartitions might be relevant
	std::set<uint32_t> crossed_subparition_indices;
	uint32_t last_subpartition_index;
	collect_crossed_subpartitions(p, vm, remaining_displacement, crossed_subparition_indices, last_subpartition_index);


	// TBD: check wall collisions
	// here we can return RAY_TRACE_HIT_WALL

	// for each SP
	for (uint32_t sp_index: crossed_subparition_indices) {


		// get cached reacting molecules for this SP
		subpartition_mask_t& sp_reactants = p.volume_molecule_reactants[sp_index][vm.species_id]; //get_sp_species_reacting_mols_cached_data(sp_index, vm, p);

		// for each molecule in this SP
		//uint32_t indices[4];
		//uint32_t pos;
		for (uint32_t vm_index: sp_reactants) {
			volume_molecule_t& colliding_vm = p.volume_molecules[vm_index];

			// we would like to compute everything that's needed just once
			float_t time;
			vec3_t position;
			// collide_mol must be inlined because many things are computed all over there
			if (collide_mol(vm, remaining_displacement, colliding_vm, time, position, world->world_constants.rx_radius_3d)) {
				reaction_t* rx = world->get_reaction(vm, colliding_vm);
				assert(rx != nullptr);
				molecule_collisions.push_back(
						molecules_collision_t(&p, vm.idx, colliding_vm.idx, rx, time, position)
				);
			}
		}
	}

	// the value is valid only when RAY_TRACE_FINISHED is returned
	new_subpartition_index = last_subpartition_index;
	new_position = vm.pos + remaining_displacement;

  return RAY_TRACE_FINISHED; // no wall was hit
}

/******************************************************************************
 *
 * collide_and_react_with_vol_mol is a helper function used in diffuse_3D to
 * handle collision of a diffusing molecule with a molecular target.
 *
 * Returns 1 if reaction does happen and 0 otherwise.
 *
 ******************************************************************************/
bool diffuse_react_event_t::collide_and_react_with_vol_mol(
		partition_t& p,
		molecules_collision_t& collision,
		vec3_t& displacement,
		float_t remaining_time_step
)	{

  volume_molecule_t& colliding_molecule = p.get_vm(collision.colliding_molecule_idx); // am
  volume_molecule_t& diffused_molecule = p.get_vm(collision.diffused_molecule_idx); // m

  // returns 1 when there are no walls at all
  //TBD: double factor = exact_disk(

  reaction_t& rx = *collision.rx;
  //  rx->prob_t is always NULL in out case update_probs(world, rx, m->t);
  // returns which reaction pathway to take
  int i = test_bimolecular(
    rx, /*scaling, 0,*/ colliding_molecule, diffused_molecule/*, world->rng*/);

  if (i < RX_LEAST_VALID_PATHWAY) {
  	return false;
  }
  else {
  	// might invalidate references
  	int j = outcome_bimolecular(p, collision, i, remaining_time_step);
  	assert(j == RX_DESTROY);
  	return true;
  }
}


/*************************************************************************
test_bimolecular
  In: the reaction we're testing
      a scaling coefficient depending on how many timesteps we've
        moved at once (1.0 means one timestep) and/or missing interaction area
      local probability factor (positive only for the reaction between two
        surface molecules, otherwise equal to zero)
      reaction partners
  Out: RX_NO_RX if no reaction occurs
       int containing which reaction pathway to take if one does occur
  Note: If this reaction does not return RX_NO_RX, then we update
        counters appropriately assuming that the reaction does take place.
*************************************************************************/
/*, double scaling - 1, double local_prob_factor = 0,*/
int diffuse_react_event_t::test_bimolecular(
		reaction_t& rx,
		volume_molecule_t& a1,
		volume_molecule_t& a2
) {
	//not true assert(a1.subpartition_index == a2.subpartition_index && "Subpartitions must be identical");

  /* rescale probabilities for the case of the reaction
     between two surface molecules */
  float_t min_noreaction_p = rx.min_noreaction_p; // local_prob_factor == 0

	/* Instead of scaling rx->cum_probs array we scale random probability */
  float_t p = rng_dbl(&(world->rng));

	if (p >= min_noreaction_p)
		return RX_NO_RX;
	else
		return 0; // we have just one pathwayy
}


int diffuse_react_event_t::outcome_bimolecular(
		partition_t& p,
		molecules_collision_t& collision,
		int path,
		float_t remaining_time_step
) {

	// might invalidate references
	int result = outcome_products_random(p, collision, path, remaining_time_step);

	if (result == RX_A_OK) {
		volume_molecule_t& reacA = p.get_vm(collision.diffused_molecule_idx);
		volume_molecule_t& reacB = p.get_vm(collision.colliding_molecule_idx);

		//assert(reacA.subpartition_index == reacB.subpartition_index && "Subpartitions must be identical");
#ifdef DEBUG_REACTIONS
		// ref. printout first destroys B then A
		DUMP_CONDITION4(
			reacB.dump(world, "", "  defunct vm:", world->current_iteration);
			reacA.dump(world, "", "  defunct vm:", world->current_iteration);
		);
#endif

		// always for now
		// we used the reactants - remove them
		p.set_molecule_as_defunct(reacA);
		p.set_molecule_as_defunct(reacB);

		return RX_DESTROY;
	}


  return result;
}

// why is this called "random"? - check if reaction occurs is in test_bimolecular
// mcell3 version returns  cross_wall ? RX_FLIP : RX_A_OK;
// ! might invalidate references
int diffuse_react_event_t::outcome_products_random(
		partition_t& p,
		molecules_collision_t& collision,
		int path,
		float_t remaining_time_step
) {
	// we can have just one product for now and no walls

	// create and place each product
	assert(collision.rx->products.size() == 1);
	volume_molecule_t vm(MOLECULE_IDX_INVALID, collision.rx->products[0].species_id, collision.position);

	volume_molecule_t& new_vm = p.add_volume_molecule(vm, world->species[vm.species_id].time_step);
	new_vm.flags =  TYPE_VOL | IN_VOLUME | ACT_DIFFUSE;

#ifdef DEBUG_REACTIONS
	DUMP_CONDITION4(
		new_vm.dump(world, "", "  created vm:", world->current_iteration);
	);
#endif

	new_molecules_to_diffuse.push_back(
			molecule_to_diffuse_t(new_vm.idx, remaining_time_step - collision.time));

	// schedule new product for diffusion - where?
	return RX_A_OK; // ?
}

void molecules_collision_t::dump(partition_t& p, const std::string ind) const {
	cout << ind << "diffused_molecule:\n";
	p.get_vm(diffused_molecule_idx).dump(ind + "  ");
	cout << ind << "colliding_molecule:\n";
	p.get_vm(colliding_molecule_idx).dump(ind + "  ");
	cout << ind << "reaction:";
	rx->dump(ind + "  ");

	cout << "time: \t\t" << time << " [float_t] \t\t\n";
	cout << "position: \t\t" << position << " [vec3_t] \t\t\n";
}

string molecules_collision_t::to_string() const {
	stringstream ss;
	ss
		//	<< "diff_idx: " << diffused_molecule_idx
			<< "coll_idx: " << colliding_molecule_idx
			<< ", time: " << time
			<< ", pos: " << position;
	return ss.str();
}

void molecules_collision_t::dump_array(partition_t& p, const std::vector<molecules_collision_t>& vec) {
	// printed in resverse - same as
	for (size_t i = 0; i < vec.size(); i++) {
		cout << "  " << "collision " << i << ": " << vec[i].to_string() << "\n";
	}
}



} /* namespace mcell */
