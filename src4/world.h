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
#include <map>


#include "partition.h"
#include "scheduler.h"
#include "species.h"
#include "reaction.h"
#include "geometry.h"

namespace mcell {


class world_t {
private:
  void init_fpu();
  void init_simulation();
  void create_defragmentation_events();

public:
  world_t();
  void init_world_constants();
  bool run_simulation();

  // -------------- parition manipulation methods --------------
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
    assert(world_constants.partition_edge_length != 0);
    assert(get_partition_index(pos) == PARTITION_INDEX_INVALID && "Partition must not exist");
    vec3_t origin =
        floor_to_multiple(pos, world_constants.partition_edge_length)
        - vec3_t(world_constants.partition_edge_length/2);
    partitions.push_back(partition_t(origin, world_constants));
    return partitions.size() - 1;
  }

  // -------------- reaction utility methods --------------

  // should be inlined
  bool can_react_vol_vol(const volume_molecule_t& a, const volume_molecule_t& b) const {
    // must not be the same molecule
    if (&a == &b) {
      return false;
    }

    // not be defunct
    if (a.is_defunct() || b.is_defunct()) {
      return false;
    }
    // is there even any reaction
    const species_t& sa = species[a.species_id];
    if (!sa.has_flag(SPECIES_FLAG_CAN_VOLVOL)) {
      return false;
    }
    const species_t& sb = species[b.species_id];
    if (!sb.has_flag(SPECIES_FLAG_CAN_VOLVOL)) {
      return false;
    }

    // is there a reaction between these two species?
    auto it_first_species = bimolecular_reactions_map.find(a.species_id);
    if (it_first_species == bimolecular_reactions_map.end()) {
      return false;
    }
    auto it_second_species = it_first_species->second.find(b.species_id);
    if (it_second_species == it_first_species->second.end()) {
      return false;
    }

    return true;
  }

  // must return result, asserts otherwise
  const reaction_t* get_reaction(const volume_molecule_t& a, const volume_molecule_t& b) const {
    const auto& it_map_for_species = bimolecular_reactions_map.find(a.species_id);
    assert(it_map_for_species != bimolecular_reactions_map.end());
    const auto& it_res = it_map_for_species->second.find(b.species_id);
    assert(it_res != it_map_for_species->second.end());
    return it_res->second;
  }

  // -------------- geometry utility methods --------------

  partition_index_t get_partition_index_for_pos(const vec3_t& pos);

  // adds vertex to a given partition and returns pair partition index and the index of the
  // vertex in that partition
  partition_vertex_index_pair_t add_geometry_vertex(const vec3_t& pos);

  // adds a new geometry object with its walls, sets unique ids for the walls and objects
  void add_geometry_object(
      const geometry_object_t& obj,
      const std::vector<wall_t>& walls,
      const std::vector<std::vector<partition_vertex_index_pair_t>>& walls_vertices
  );


  void dump();

  // -------------- world data --------------
  std::vector<partition_t> partitions;

  scheduler_t scheduler;

  std::vector<species_t> species; // owner

  std::vector<reaction_t> reactions; // we might need faster searching or reference from species to reactions here but let's keep it simple for now

  // FIXME: there might be multiple reactions for 1 or 2 reactants (multiple pathways)
  unimolecular_reactions_map_t unimolecular_reactions_map; // created from reactions in init_simulation
  bimolecular_reactions_map_t bimolecular_reactions_map; // created from reactions in init_simulation

  uint64_t current_iteration;
  uint64_t iterations; // number of iterations to simulate


  uint32_t seed_seq;

  world_constants_t world_constants;

  // single state for the random number generator
  rng_state rng;

  // in case when there would be many copies of a string, this constant pool can be used
  const char* add_const_string_to_pool(const std::string str) {
    return const_string_pool.insert(str).first->c_str();
  }
private:
  std::set<std::string> const_string_pool;

  wall_id_t next_wall_id;
  geometry_object_id_t next_geometry_object_id;
};

} // namespace mcell

#endif // SRC4_WORLD_H_

