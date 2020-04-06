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

#include "../libs/bng/rxn_container.h"
#include "bng/bng.h"

#include "partition.h"
#include "scheduler.h"
#include "geometry.h"
#include "callback_info.h"
#include "count_buffer.h"
#include "counted_volumes_util.h"


namespace MCell {


class World {
public:
  World();
  void init_simulation();
  void run_simulation(const bool dump_initial_state = false);
  void run_n_iterations(
      const uint64_t num_iterations,
      const uint64_t output_frequency,
      const bool terminate_last_iteration_after_viz_output = false // needed for exact match with MCell3, must false when used from pymcell
  );
  void end_simulation();

  // -------------- diverse getters -----------------------------
  const SimulationConfig& get_config() {
    return config;
  }

  uint64_t get_current_iteration() const {
    return stats.get_current_iteration();
  }

  // -------------- partition manipulation methods --------------
  partition_id_t get_partition_index(const Vec3& pos) {
    // for now a slow approach, later some hashing/memoization might be needed
    for (partition_id_t i = 0; i < partitions.size(); i++) {
      if (partitions[i].in_this_partition(pos)) {
        return i;
      }
    }
    return PARTITION_ID_INVALID;
  }

  partition_id_t get_or_add_partition_index(const Vec3& pos) {

    partition_id_t res = get_partition_index(pos);
    // not found - add a new partition
    if (res == PARTITION_ID_INVALID) {
      res = add_partition(pos);
    }

    return res;
  }

  // add a partition in a predefined 'lattice' that contains point pos
  partition_id_t add_partition(const Vec3& pos) {
    assert(config.partition_edge_length != 0);
    assert(get_partition_index(pos) == PARTITION_ID_INVALID && "Partition must not exist");

    Vec3 origin =
        floor_to_multiple(pos, config.partition_edge_length)
        - Vec3(config.partition_edge_length/2);

    partitions.push_back(Partition(partitions.size(), origin, config, bng_engine, stats));
    return partitions.size() - 1;
  }

  Partition& get_partition(partition_id_t i) {
    assert(i < partitions.size());
    return partitions[i];
  }

  const Partition& get_partition(partition_id_t i) const {
    assert(i < partitions.size());
    return partitions[i];
  }

  PartitionVector& get_partitions() {
      return partitions;
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

  counted_volume_id_t get_next_counted_volume_id() {
    counted_volume_id_t res = next_counted_volume_id;
    next_counted_volume_id++;
    return res;
  }

  void dump();

  // -------------- callback registration --------------
  void register_wall_hit_callback_internal(wall_hit_callback_func func, void* clientdata_) {
    wall_hit_callback = func;
    wall_hit_callback_clientdata = clientdata_;
  }

  wall_hit_callback_func get_wall_hit_callback() {
    return wall_hit_callback;
  }

  // ---------------------- other ----------------------
  BNG::SpeciesContainer& get_all_species() { return bng_engine.get_all_species(); }
  const BNG::SpeciesContainer& get_all_species() const { return bng_engine.get_all_species(); }

  BNG::RxnContainer& get_all_rxns() { return bng_engine.get_all_rxns(); }
  const BNG::RxnContainer& get_all_rxns() const { return bng_engine.get_all_rxns(); }

  count_buffer_id_t create_count_buffer(const std::string filename, const size_t buffer_size) {
    count_buffer_id_t id = count_buffers.size();
    count_buffers.push_back(CountBuffer(filename, buffer_size));
    return id;
  }

  CountBuffer& get_count_buffer(const count_buffer_id_t id) {
    assert(id < count_buffers.size());
    return count_buffers[id];
  }


private:
  void init_fpu();
  void create_defragmentation_events();
  void init_counted_volumes();

public:
  // single instance for the whole mcell simulator,
  // used as constants during simulation
  SimulationConfig config;

  BNG::BNGEngine bng_engine;

  SimulationStats stats;

  Scheduler scheduler;

  uint64_t iterations; // number of iterations to simulate - move to Sim config
  uint seed_seq; // initial seed passed to mcell as argument

  rng_state rng; // single state for the random number generator

private:
  PartitionVector partitions;

  CountBufferVector count_buffers;

  // global ID counters
  wall_id_t next_wall_id;
  geometry_object_id_t next_geometry_object_id;
  counted_volume_id_t next_counted_volume_id;

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

public:
  // NOTE: only a temporary solution of callbacks for now

  // callbacks
  wall_hit_callback_func wall_hit_callback;
  // clientdata hold information on what Python function we should call
  void* wall_hit_callback_clientdata;
};

} // namespace mcell

#endif // SRC4_WORLD_H_

