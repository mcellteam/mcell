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

#include <time.h>
#include <sys/time.h> // Linux include
#include <sys/resource.h> // Linux include

#include <vector>
#include <string>
#include <set>
#include <map>

#include "world_constants.h"
#include "partition.h"
#include "scheduler.h"
#include "species.h"
#include "reaction.h"
#include "geometry.h"
#include "callback_info.h"

namespace MCell {


class World {
private:
  void init_fpu();
  void init_simulation();
  void create_defragmentation_events();

public:
  World();
  void init_world_constants();
  void run_simulation(const bool dump_initial_state = false);
  void run_n_iterations(const uint64_t num_iterations, const uint64_t output_frequency);
  void end_simulation();

  // -------------- partition manipulation methods --------------
  partition_index_t get_partition_index(const vec3_t& pos) {
    // for now a slow approach, later some hashing/memoization might be needed
    for (partition_index_t i = 0; i < partitions.size(); i++) {
      if (partitions[i].in_this_partition(pos)) {
        return i;
      }
    }
    return PARTITION_INDEX_INVALID;
  }

  partition_index_t get_or_add_partition_index(const vec3_t& pos) {

    partition_index_t res = get_partition_index(pos);
    // not found - add a new partition
    if (res == PARTITION_INDEX_INVALID) {
      res = add_partition(pos);
    }

    return res;
  }

  // add a partition in a predefined 'lattice' that contains point pos
  partition_index_t add_partition(const vec3_t& pos) {
    assert(world_constants.partition_edge_length != 0);
    assert(get_partition_index(pos) == PARTITION_INDEX_INVALID && "Partition must not exist");

    vec3_t origin =
        floor_to_multiple(pos, world_constants.partition_edge_length)
        - vec3_t(world_constants.partition_edge_length/2);

    partitions.push_back(Partition(origin, world_constants, simulation_stats));
    return partitions.size() - 1;
  }

  Partition& get_partition(partition_index_t i) {
    assert(i < partitions.size());
    return partitions[i];
  }

  std::vector<Partition>& get_partitions() {
      return partitions;
  }

  // -------------- reaction utility methods --------------

  // should be inlined
  bool can_react_vol_vol(const Molecule& a, const Molecule& b) const {
    // must not be the same molecule
    if (&a == &b) {
      return false;
    }

    // not be defunct
    if (a.is_defunct() || b.is_defunct()) {
      return false;
    }
    // is there even any reaction
    const Species& sa = species[a.species_id];
    if (!sa.has_flag(SPECIES_FLAG_CAN_VOLVOL)) {
      return false;
    }
    const Species& sb = species[b.species_id];
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
  const Reaction* get_reaction(const Molecule& a, const Molecule& b) const {
    const auto& it_map_for_species = bimolecular_reactions_map.find(a.species_id);
    assert(it_map_for_species != bimolecular_reactions_map.end());
    const auto& it_res = it_map_for_species->second.find(b.species_id);
    assert(it_res != it_map_for_species->second.end());
    return it_res->second;
  }

  const Species& get_species(const species_id_t species_id) const {
    assert(species_id < species.size());
    return species[species_id];
  }

  const std::vector<Species>& get_species() const {
    return species;
  }

  void add_species(const Species& new_species) {
    assert(new_species.species_id == species.size());
    species.push_back(new_species);
  }
  // -------------- object id counters --------------
  wall_id_t get_next_wall_id() {
    wall_id_t res = next_wall_id;
    next_wall_id++;
    return res;
  }

  geometry_object_id_t get_next_geometry_object_id() {
    geometry_object_id_t res = next_geometry_object_id;
    next_geometry_object_id++;
    return res;
  }

  void dump();

  // -------------- callback registration --------------
  void register_wall_hit_callback_internal(wall_hit_callback_func func, void* clientdata_) {
    wall_hit_callback = func;
    wall_hit_callback_clientdata = clientdata_;
  }

  // clientdata hold information on what Python function we should call
  void* wall_hit_callback_clientdata;

  wall_hit_callback_func get_wall_hit_callback() {
    return wall_hit_callback;
  }

private:
  std::vector<Partition> partitions;
  std::vector<Species> species;

public:
  Scheduler scheduler;

  std::vector<Reaction> reactions; // we might need faster searching or reference from species to reactions here but let's keep it simple for now

  // TODO_PATHWAYS: there might be multiple reactions for 1 or 2 reactants (multiple pathways)
  UnimolecularReactionsMap unimolecular_reactions_map; // created from reactions in init_simulation
  BimolecularReactionsMap bimolecular_reactions_map; // created from reactions in init_simulation

  uint64_t current_iteration;
  uint64_t iterations; // number of iterations to simulate

  uint seed_seq; // initial seed passed to mcell as argument

  WorldConstants world_constants;
  SimulationStats simulation_stats;

  rng_state rng; // single state for the random number generator

  // in case when there would be many copies of a string, this constant pool can be used
  const char* add_const_string_to_pool(const std::string str) {
    return const_string_pool.insert(str).first->c_str();
  }

private:
  std::set<std::string> const_string_pool;

  // global ID counters
  wall_id_t next_wall_id;
  geometry_object_id_t next_geometry_object_id;


  // used by run_n_iterations to know whether the simulation was
  // already initialized
  bool simulation_initialized;

  // used to know whether we already reported final simulation stats
  bool simulation_ended;

  // several variables to report simulation time
  timeval previous_progress_report_time;
  rusage sim_start_time;

  // and to nicely report simulation progress
  uint64_t previous_iteration;

  // callbacks
  wall_hit_callback_func wall_hit_callback;
};

} // namespace mcell

#endif // SRC4_WORLD_H_

