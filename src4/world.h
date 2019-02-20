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

#ifndef SRC4_WORLD_H_
#define SRC4_WORLD_H_

#include <vector>
#include <string>
#include <set>

#include "partition.h"
#include "scheduler.h"
#include "species.h"

namespace mcell {


class world_t {
private:
	void init();
public:
	bool run_simulation();

	uint32_t get_partition_index(const vec3_t& pos) {
		// for now a slow approach, later some hashing/memoization might be needed
		for (uint32_t i = 0; i < partitions.size(); i++) {
			if (partitions[i].in_this_partition(pos)) {
				return i;
			}
		}
		return PARTITION_INDEX_INVALID;
	}

	uint32_t get_or_add_partition_index(const vec3_t& pos) {
		uint32_t res;
		res = get_partition_index(pos);
		if (res == PARTITION_INDEX_INVALID) {
			res = add_partition(pos);
		}
		// not found - add a new partition
		return res;
	}

	// add a partition in a predefined 'lattice' that contains point pos
	uint32_t add_partition(const vec3_t& pos) {
		// TODO: some check on validity of pos?
		assert(get_partition_index(pos) == PARTITION_INDEX_INVALID && "Parition must not exist");

		vec3_t origin = floor_to_multiple(pos, partition_edge_length) - vec3_t(partition_edge_length/2);

		partitions.push_back(partition_t(origin, partition_edge_length));
		partition_t& new_partition = partitions.back();

		return partitions.size() - 1;
	}

  std::vector<partition_t> partitions;

  scheduler_t scheduler;

  std::vector<species_t> species;

  float_t time_unit;
  uint64_t iterations; // number of iterations to simulate

  float_t length_unit;

  uint32_t seed_seq;

  float_t partition_edge_length;

  // single state for the random number generator
  rng_state rng;

  // in case when there would be many copies of a string, this constant pool can be used
  const char* add_const_string_to_pool(const std::string str) {
  	return const_string_pool.insert(str).first->c_str();
  }
private:
  std::set<std::string> const_string_pool;
};

} /* namespace mcell4 */

#endif /* SRC4_WORLD_H_ */
