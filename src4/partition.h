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

#include "defines.h"
#include "molecule.h"

namespace mcell {

// do not call set() directly, maybe there is a way how to forbid this
// TODO: use a different representation, maybe a set for now

class subpartition_mask_t: public std::unordered_set<uint32_t> {
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
};

#if 0
class subpartition_mask_t: public std::bitset<MAX_MOLECULES_PER_PARTITION> {
public:
	// only adding a check to the inherited set method
	void set_contains_molecule(uint32_t index, bool value = true) {
		assert(index < size());
		set(index, value);
		indexing_initialized = false;
	}

	// note: maybe change into iterators
	uint32_t get_first_index() {
		indexing_initialized = true;

		for (uint32_t i = 0; i < size(); i++) {
			if (this->operator [](i)) {
				last_returned_index = i;
				return i;
			}
		}
		last_returned_index = MOLECULE_INDEX_INVALID;
		return MOLECULE_INDEX_INVALID;
	}

	uint32_t get_next_index() {
		assert(indexing_initialized);

		for (uint32_t i = last_returned_index + 1; i < size(); i++) {
			if (this->operator [](i)) {
				last_returned_index = i;
				return i;
			}
		}
		last_returned_index = MOLECULE_INDEX_INVALID;
		return MOLECULE_INDEX_INVALID;
	}

private:
	bool indexing_initialized;
	uint32_t last_returned_index;
};
#endif

class partition_t {
public:
	partition_t(const vec3_t origin_, const world_constants_t& world_constants_)
		: origin_corner(origin_),
			world_constants(world_constants_) {
		opposite_corner = origin_corner + world_constants.partition_edge_length;
		// preaallocate volume_molecules arrays and also volume_molecule_indices_per_time_step
		volume_molecules_subpartition_masks.resize(powu(world_constants.subpartitions_per_partition_dimension, 3));
	}


	void dump();

	uint32_t get_molecule_index(const volume_molecule_t& m) {
		// simply use pointer arithmetic to compute the molecule's index
		int res = &m - &volume_molecules[0];
		assert(res >= 0);
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
				pair_time_step_volume_molecules_t(time_step, std::vector< molecule_index_t >()));
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
		res = relative_position / world_constants.subpartition_edge_length;
	}

	uint32_t get_subpartition_index_from_3d_indices(const ivec3_t& indices) const {
		uint32_t dim = world_constants.subpartitions_per_partition_dimension;
		// volume_molecules_subpartition_masks is a flattened cube of dimension dim
		// example: dim: 5x5x5,  (1, 2, 3) -> 1 + 2*5 + 3*5*5 = 86
		return indices.x + indices.y * dim + indices.z * powu(dim, 2);
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

	void change_molecule_subpartition(volume_molecule_t& m, const uint32_t new_subpartition_index) {
		assert(m.subpartition_index < volume_molecules_subpartition_masks.size());
		assert(new_subpartition_index < volume_molecules_subpartition_masks.size());
		if (m.subpartition_index == new_subpartition_index) {
			return; // nothing to do
		}
		uint32_t molecule_index = get_molecule_index(m);
#ifdef DEBUG_PARTITION
		std::cout << "Molecule " << molecule_index << " changed subpartition from "
				<<  m.subpartition_index << " to " << new_subpartition_index << ".\n";
#endif
		// clear old position and set a new one
		volume_molecules_subpartition_masks[m.subpartition_index].set_contains_molecule(molecule_index, false);
		volume_molecules_subpartition_masks[new_subpartition_index].set_contains_molecule(molecule_index, true);
		m.subpartition_index = new_subpartition_index;
	}

	void add_volume_molecule(const volume_molecule_t& m, const uint32_t time_step_index) {
		uint32_t molecule_index = volume_molecules.size();
		// and its index to the list sorted by time step
		// this is an array that changes only when molecule leaves this partition
		assert(time_step_index <= volume_molecule_indices_per_time_step.size());
		volume_molecule_indices_per_time_step[time_step_index].second.push_back(molecule_index);

		// compute subpartition and set it, this might change often, therefore we are using bitfield that
		// can be easily changed
		uint32_t subpartition_index = get_subpartition_index(m.pos);
		assert(subpartition_index < volume_molecules_subpartition_masks.size());
		volume_molecules_subpartition_masks[subpartition_index].set_contains_molecule(molecule_index, true);

		// and finally store (copy) the new molecule and set its subpartition index
		volume_molecules.push_back(m);
		volume_molecules.back().subpartition_index = subpartition_index;
	}

	// left, bottom, closest (lowest z) point of the partition
  vec3_t origin_corner;
  vec3_t opposite_corner;

  // vector containing all volume molecules in this partition
  std::vector< /* molecule index*/ volume_molecule_t> volume_molecules;

  // arrays of indices to the volume_molecules array where each array corresponds to a given time step
  typedef std::pair< float_t, std::vector< molecule_index_t > > pair_time_step_volume_molecules_t;
  // indexed by diffusion time step index
  std::vector< pair_time_step_volume_molecules_t > volume_molecule_indices_per_time_step;

  // indexed by subpartition index, size is world->subpartition_edge_length^3
  std::vector < subpartition_mask_t > volume_molecules_subpartition_masks;

  //TBD: std::vector< /* surface molecule index */ surface_molecule> surface_molecules;
  //TBD: std::vector< /* subpartition index */ subpartition_mask > surface_molecules_subpatition_masks;

  const world_constants_t& world_constants;
};

} // namespace mcell

#endif /* SRC4_PARTITION_H_ */
