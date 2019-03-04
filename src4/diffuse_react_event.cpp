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
	cache_sp_species_reacting_mols.clear();

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

	// Possible optimization:
	// GPU version - can make sets of possible collisions for a given species before
	// the whole diffuse is run
	// ! CPU version - compute sets of possible colilisions for a given species and "cache" it

	for (molecule_index_t i: indices) {
		volume_molecule_t& vm = p.volume_molecules[i];
		if (vm.is_defunct())
			continue;
		species_t& spec = world->species[vm.species_id];


/*	  double steps = 1.0;
	  double t_steps = 1.0;
	  double rate_factor = 1.0;
	  double r_rate_factor = 1.0;
*/

		// 2) diffuse each molecule - get information on position change
		// TBD: reflections
		vec3_t displacement;
		compute_displacement(spec, displacement);

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

#ifdef DEBUG_COLLISIONS
			molecules_collision_t::dump_array(molecule_collisions);
#endif

		} while (state != RAY_TRACE_FINISHED);

		// 4) evaluate and possible execute reactions
		// are they sorted by time?
		for (molecules_collision_t& collision: molecule_collisions) {
			// reactions with time 0 are ignored


			assert(collision.time > 0 && collision.time <= 1);

			// mol-mol collision
			if (collision.time < EPS) {
					continue;
			}

			/*// ??
	    struct vector3* loc_certain = NULL;
	    struct collision *tentative = shead2; // pointer to the first collision

		  //double t_steps = 1.0;
		  //double r_rate_factor = 1.0;

			collide_and_react_with_vol_mol(
					collision,
			*/

		}

		// TODO
		/*
		 *       if ((smash->what & COLLIDE_VOL) != 0) {
        if (smash->t < EPS_C) {
          continue;
        }
        if (collide_and_react_with_vol_mol(world, smash, vm, &tentative,
          &displacement, loc_certain, t_steps, r_rate_factor) == 1) {
          FREE_COLLISION_LISTS();
          return NULL;
        } else {
          continue;
        }
      } else if ((smash->what & COLLIDE_WALL) != 0) {
		 */

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

// TODO: optimize
subpartition_mask_t& diffuse_react_event_t::get_sp_species_reacting_mols_cached_data(
		uint32_t sp_index, volume_molecule_t& vm, partition_t& p) {

	auto it_per_species = cache_sp_species_reacting_mols.find(sp_index);
	if (it_per_species != cache_sp_species_reacting_mols.end()) {
		auto it_per_sp_mask = it_per_species->second.find(vm.species_id);
		if (it_per_sp_mask != it_per_species->second.end()) {
			return it_per_sp_mask->second;
		}

	}

	// not found
	std::pair<species_id_t, subpartition_mask_t> new_cache_item(vm.species_id, subpartition_mask_t());
	auto it_cached = cache_sp_species_reacting_mols[sp_index].insert(new_cache_item);
	subpartition_mask_t& mask_to_init = it_cached.first->second;

	subpartition_mask_t& full_sp_mask = p.volume_molecules_subpartition_masks[sp_index];

	for (uint32_t vm_index: full_sp_mask) {
				volume_molecule_t& colliding_vm = p.volume_molecules[vm_index];

				// can we react?
				if (world->can_react_vol_vol(vm, colliding_vm)) {
					mask_to_init.insert(vm_index);
				}
	}
	return mask_to_init;
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


		// get cached reacting molecules for this SP
		subpartition_mask_t& sp_cached_mask = get_sp_species_reacting_mols_cached_data(sp_index, vm, p);

		// for each molecule in this SP
		for (uint32_t vm_index: sp_cached_mask) {
			volume_molecule_t& colliding_vm = p.volume_molecules[vm_index];

			// can we react? - already pre-cached
			//if (!world->can_react_vol_vol(vm, colliding_vm)) {
			//	continue;
			//}

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

/******************************************************************************
 *
 * collide_and_react_with_vol_mol is a helper function used in diffuse_3D to
 * handle collision of a diffusing molecule with a molecular target.
 *
 * Returns 1 if reaction does happen and 0 otherwise.
 *
 ******************************************************************************/
#if 0
static int collide_and_react_with_vol_mol(
	struct volume* world,
  struct collision* smash,
	struct volume_molecule* m,
	struct collision** tentative,
	struct vector3* displacement,
	struct vector3* loc_certain, double  t_steps, double r_rate_factor) {

  struct abstract_molecule* am = (struct abstract_molecule *)smash->target;
  double factor = exact_disk(
      world, &(smash->loc), displacement, world->rx_radius_3d, m->subvol,
      m, (struct volume_molecule *)am, world->use_expanded_list,
      world->x_fineparts, world->y_fineparts, world->z_fineparts);

  if (factor < 0) { /* Probably hit a wall, might have run out of memory */
    return 0; /* Reaction blocked by a wall */
  }

  double scaling = factor * r_rate_factor;
  struct rxn* rx = smash->intermediate;
  if ((rx != NULL) && (rx->prob_t != NULL)) {
    update_probs(world, rx, m->t);
  }

  struct species *spec = m->properties;
  struct periodic_image *periodic_box = m->periodic_box;
  int i = test_bimolecular(
    rx, scaling, 0, am, (struct abstract_molecule *)m, world->rng);

  if (i < RX_LEAST_VALID_PATHWAY) {
    return 0;
  }

  int j = outcome_bimolecular(world, rx, i, (struct abstract_molecule *)m, am,
    0, 0, m->t + t_steps * smash->t, &(smash->loc), loc_certain);

  if (j != RX_DESTROY) {
    return 0;
  } else {
    /* Count the hits up until we were destroyed */
    struct collision* ttv = *tentative;
    for (; ttv != NULL && ttv->t <= smash->t; ttv = ttv->next) {
      if (!(ttv->what & COLLIDE_WALL)) {
        continue;
      }
      if (m->properties == NULL) {
        continue;
      }
      if (!(m->properties->flags & ((struct wall *)ttv->target)->flags &
            COUNT_SOME_MASK)) {
        continue;
      }

      count_region_update(world, m, spec, m->id, periodic_box,
        ((struct wall *)ttv->target)->counting_regions,
        ((ttv->what & COLLIDE_MASK) == COLLIDE_FRONT) ? 1 : -1, 0, &(ttv->loc), ttv->t);
      if (ttv == smash) {
        break;
      }
    }
    *tentative = ttv;
  }
  return 1;
}
#endif


void molecules_collision_t::dump(const std::string ind) const {
	cout << ind << "diffused_molecule:";
	diffused_molecule.dump(ind + "  ");
	cout << ind << "colliding_molecule:";
	colliding_molecule.dump(ind + "  ");

	cout << "time: \t\t" << time << " [float_t] \t\t\n";
	cout << "position: \t\t" << position << " [vec3_t] \t\t\n";
}

void molecules_collision_t::dump_array(const std::vector<molecules_collision_t>& vec) {
	cout << "Collision array: " << (vec.empty() ? "EMPTY" : "") << "\n";

	for (size_t i = 0; i < vec.size(); i++) {
		cout << i << ":\n";
		vec[i].dump("  ");
	}
}



} /* namespace mcell */
