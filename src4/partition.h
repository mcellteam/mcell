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

#include "defines.h"
#include "molecule.h"
#include "molecule.h"

namespace mcell {

class partition_t {
public:
	partition_t(const vec3_t origin_, const float_t& partition_edge_length)
		: origin_corner(origin_) {
		opposite_corner = origin_corner + partition_edge_length;
		// TODO: preaallocate volume_molecules arrays and
		// also volume_molecule_indices_per_time_step
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

	bool in_this_partition(const vec3_t& pos) {
		return glm::all(glm::greaterThanEqual(pos, origin_corner))
			&& glm::all(glm::lessThan(pos, opposite_corner));
	}


	void add_volume_molecule(const volume_molecule_t& m, const uint32_t time_step_index) {
		// store the new molecule
		volume_molecules.push_back(m);
		// and also its index
		assert(time_step_index <= volume_molecule_indices_per_time_step.size());
		volume_molecule_indices_per_time_step[time_step_index].second.push_back(volume_molecules.size() - 1);
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

  // TBD
  //std::vector< /* subpartition index */
  //std::vector < /* diffusion time step index */ subpartition_mask_t > > volume_molecules_subpartition_masks;

  //TBD: std::vector< /* surface molecule index */ surface_molecule> surface_molecules;
  //TBD: std::vector< /* subpartition index */ subpartition_mask > surface_molecules_subpatition_masks;

};

} // namespace mcell

#endif /* SRC4_PARTITION_H_ */
