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

#ifndef SRC4_PARTITION_H_
#define SRC4_PARTITION_H_

#include <unordered_set>
#include <set>

#include "defines.h"
#include "molecule.h"

namespace mcell {

// do not call set() directly, maybe there is a way how to forbid this
// TODO: use a different representation, maybe a set for now

// FIXME: use molecule_idx_t
class subpartition_mask_t: public std::set<uint32_t>  // set is faster than unordered_set
		{
public:
	void set_contains_molecule(uint32_t index, bool value = true) {
		if (value) {
			assert(count(index) == 0);
			insert(index);
		}
		else {
			assert(count(index) == 1);
			erase(index);
		}
	}

	void dump();
};

class partition_t {
public:
	partition_t(const vec3_t origin_, const world_constants_t& world_constants_)
		: origin_corner(origin_),
			next_molecule_idx(0),
			world_constants(world_constants_) {

		opposite_corner = origin_corner + world_constants.partition_edge_length;
		// preaallocate volume_molecules arrays and also volume_molecule_indices_per_time_step
		uint32_t num_subparts = powu(world_constants.subpartitions_per_partition_dimension, 3);
		volume_molecule_reactants.resize(num_subparts);

		size_t num_species = world_constants.bimolecular_reactions_map->size();
		for (auto& reactants : volume_molecule_reactants) {
			reactants.resize(num_species);
		}

	}


	void dump();

	volume_molecule_t& get_vm(const molecule_idx_t idx) {
		assert(idx != MOLECULE_IDX_INVALID && idx < volume_molecules.size());
		return volume_molecules[idx];
	}


	molecule_idx_t get_molecule_index(const volume_molecule_t& m) {
		// simply use pointer arithmetic to compute the molecule's index
		molecule_idx_t res = m.idx;
		assert(res != MOLECULE_IDX_INVALID);
		return res;
	}

	uint32_t get_molecule_list_index_for_time_step(const float_t time_step) {
		// TODO: memoization
		for (uint32_t i = 0; i < volume_molecule_indices_per_time_step.size(); i++) {
			if (volume_molecule_indices_per_time_step[i].first == time_step) {
				return i;
			}
		}
		return PARTITION_INDEX_INVALID;
	}

	uint32_t get_or_add_molecule_list_index_for_time_step(const float_t time_step) {
		uint32_t res;
		res = get_molecule_list_index_for_time_step(time_step);
		if (res == PARTITION_INDEX_INVALID) {
			volume_molecule_indices_per_time_step.push_back(
				pair_time_step_volume_molecules_t(time_step, std::vector< molecule_idx_t >()));
			res = volume_molecule_indices_per_time_step.size() - 1;
		}
		return res;
	}

	bool in_this_partition(const vec3_t& pos) const {
		return glm::all(glm::greaterThanEqual(pos, origin_corner))
			&& glm::all(glm::lessThan(pos, opposite_corner));
	}

	void get_subpartition_3d_indices(const vec3_t& pos, ivec3_t& res) const {
		assert(in_this_partition(pos));
		vec3_t relative_position = pos - origin_corner;
		//res = relative_position / world_constants.subpartition_edge_length;
		res = relative_position * world_constants.subpartition_edge_length_rcp;
	}

	uint32_t get_subpartition_index_from_3d_indices(const ivec3_t& indices) const {
		uint32_t dim = world_constants.subpartitions_per_partition_dimension;
		// example: dim: 5x5x5,  (1, 2, 3) -> 1 + 2*5 + 3*5*5 = 86
		return indices.x + indices.y * dim + indices.z * powu(dim, 2);
	}

	uint32_t get_subpartition_index_from_3d_indices(const int x, const int y, const int z) const {
		uint32_t dim = world_constants.subpartitions_per_partition_dimension;
		// example: dim: 5x5x5,  (1, 2, 3) -> 1 + 2*5 + 3*5*5 = 86

		// check for calls from collect_neigboring_subparitions, can occur, must be fixed
		assert(x >= 0);
		assert(y >= 0);
		assert(z >= 0);

		return x + y * dim + z * powu(dim, 2);
	}

	void get_subpartition_3d_indices_from_index(const uint32_t index, ivec3_t& indices) const {
		uint32_t dim = world_constants.subpartitions_per_partition_dimension;
		// example: dim: 5x5x5,  86 -> (86%5, (86/5)%5, (86/(5*5))%5) = (1, 2, 3)
		indices.x = index % dim;
		indices.y = (index / dim) % dim;
		indices.z = (index / powu(dim, 2)) % dim;
	}

	uint32_t get_subpartition_index(const vec3_t& pos) {
		ivec3_t indices;
		get_subpartition_3d_indices(pos, indices);
		return get_subpartition_index_from_3d_indices(indices);
	}

public:
	void change_reactants_map(volume_molecule_t& vm, const uint32_t new_subpartition_index, bool adding, bool removing) {
		assert(world_constants.bimolecular_reactions_map->find(vm.species_id) != world_constants.bimolecular_reactions_map->end());

		// these are all the sets of indices of reactants for this particular subpartition
		species_reactants_map_t& subpart_reactants_orig_sp = volume_molecule_reactants[vm.subpartition_index];
		species_reactants_map_t& subpart_reactants_new_sp = volume_molecule_reactants[new_subpartition_index];

		// and these are indices of possible reactants with our reactant_species_id
		const species_reaction_map_t& reactions_map = world_constants.bimolecular_reactions_map->find(vm.species_id)->second;

		// we need to set/clear flag that says that second_reactant_info.first can react with reactant_species_id
		for (const auto& second_reactant_info : reactions_map) {
			if (removing) {
				subpart_reactants_orig_sp[second_reactant_info.first].set_contains_molecule(vm.idx, false);
			}
			if (adding) {
				subpart_reactants_new_sp[second_reactant_info.first].set_contains_molecule(vm.idx, true);
			}
		}
	}

	void change_molecule_subpartition(volume_molecule_t& vm, const uint32_t new_subpartition_index) {
		assert(vm.subpartition_index < volume_molecule_reactants.size());
		assert(new_subpartition_index < volume_molecule_reactants.size());
		if (vm.subpartition_index == new_subpartition_index) {
			return; // nothing to do
		}
#ifdef DEBUG_SUBPARTITIONS
		std::cout << "Molecule " << molecule_idx << " changed subpartition from "
				<<  vm.subpartition_index << " to " << new_subpartition_index << ".\n";
#endif

		change_reactants_map(vm, new_subpartition_index, true, true);
		vm.subpartition_index = new_subpartition_index;
	}

	// version that computes the right time_step_index each time it is called
	volume_molecule_t& add_volume_molecule(const volume_molecule_t& vm_copy, const float_t time_step) {
		uint32_t time_step_index = get_or_add_molecule_list_index_for_time_step(time_step);
		return add_volume_molecule_with_time_step_index(vm_copy, time_step_index);
	}


	volume_molecule_t& add_volume_molecule_with_time_step_index(volume_molecule_t vm_copy, const uint32_t time_step_index) {
		uint32_t molecule_index = next_molecule_idx;
		next_molecule_idx++;
		// and its index to the list sorted by time step
		// this is an array that changes only when molecule leaves this partition
		assert(time_step_index <= volume_molecule_indices_per_time_step.size());
		volume_molecule_indices_per_time_step[time_step_index].second.push_back(molecule_index);

		volume_molecules.push_back(vm_copy);
		volume_molecule_t& new_vm = volume_molecules.back();

		new_vm.idx = molecule_index;
		new_vm.subpartition_index = get_subpartition_index(vm_copy.pos);

		change_reactants_map(new_vm, new_vm.subpartition_index, true, false);

		return new_vm;
	}

	void set_molecule_as_defunct(volume_molecule_t& vm) {
		// set that this molecule does not exist anymore
		vm.set_is_defunct();

		// we will keep it in diffusion arrays (volume_molecule_indices_per_time_step)
		// but we should remove it from subpartition mask (although there are check in the collision detection code for that)
		uint32_t molecule_index = get_molecule_index(vm);

		change_reactants_map(vm, 0/*unused*/, false, true);
	}

	// left, bottom, closest (lowest z) point of the partition
  vec3_t origin_corner;
  vec3_t opposite_corner;

  // vector containing all volume molecules in this partition
  std::vector< /* molecule idx*/ volume_molecule_t> volume_molecules;
  molecule_idx_t next_molecule_idx;

  // arrays of indices to the volume_molecules array where each array corresponds to a given time step
  typedef std::pair< float_t, std::vector< molecule_idx_t > > pair_time_step_volume_molecules_t;
  // indexed by diffusion time step index
  std::vector< pair_time_step_volume_molecules_t > volume_molecule_indices_per_time_step; // TODO: rename so that the name has something to do with diffusion? diffusion list?

  // indexed by subpartition index, size is world->subpartition_edge_length^3
  //std::vector < subpartition_mask_t > volume_molecules_subpartition_masks;

  typedef std::vector < /* indexed with species_id_t, we might need some hash later*/ subpartition_mask_t > species_reactants_map_t;
  std::vector < /*subpartition index*/species_reactants_map_t > volume_molecule_reactants;

  //TBD: std::vector< /* surface molecule index */ surface_molecule> surface_molecules;
  //TBD: std::vector< /* subpartition index */ subpartition_mask > surface_molecules_subpatition_masks;

  const world_constants_t& world_constants;
};

} // namespace mcell

#endif /* SRC4_PARTITION_H_ */
