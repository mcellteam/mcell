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
#include <functional>

#include <vector>
#include <string>
#include <set>
#include <map>
#include <iostream>

#include "../libs/bng/rxn_container.h"
#include "bng/bng.h"

#include "partition.h"
#include "scheduler.h"
#include "geometry.h"
#include "count_buffer.h"
#include "counted_volumes_util.h"
#include "legacy_callback_info.h"
#include "species_flags_analyzer.h"
#include "memory_limit_checker.h"

#include "logging.h"


namespace Json {
class Value;
}

namespace MCell {

namespace API {
class Model;
class Callbacks;
}

class World {

public:
  World(API::Callbacks& callbacks_);
  ~World();
  // MCell MDL mode
  void init_and_run_simulation(const bool dump_initial_state = false, const bool dump_with_geometry = false);

  void init_simulation(const float_t start_time);

  // MCell Python mode, init_simulation must be called first
  // returns number of executed iterations
  uint64_t run_n_iterations(
      const uint64_t num_iterations,
      const bool terminate_last_iteration_after_viz_output = false // used when ending simulation
  );
  void end_simulation(const bool print_final_report = true);

  // used by converters
  void create_initial_surface_region_release_event();

  // -------------- diverse getters -----------------------------
  const SimulationConfig& get_config() {
    return config;
  }

  uint64_t get_current_iteration() const {
    return stats.get_current_iteration();
  }

  API::Callbacks& get_callbacks() {
    return callbacks;
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

  // add a partition in a predefined 'lattice' that contains point pos as its llf point
  // size is given by config
  partition_id_t add_partition(const Vec3& partition_llf) {
    assert(config.partition_edge_length != 0);
    assert(get_partition_index(partition_llf) == PARTITION_ID_INVALID && "Partition must not exist");
    partitions.push_back(Partition(partitions.size(), partition_llf, config, bng_engine, stats, species_flags_analyzer));
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

  region_id_t get_next_region_id() {
    region_id_t res = next_region_id;
    next_region_id++;
    return res;
  }

  geometry_object_id_t get_next_geometry_object_id() {
    geometry_object_id_t res = next_geometry_object_id;
    next_geometry_object_id++;
    return res;
  }

  void print_periodic_stats() const;

  void dump(const bool with_geometry = false);

  // returns empty string if everything went well, nonempty string with error message
  std::string export_to_bngl(const std::string& file_name) const;

  // the export to directory is usually called periodically and the output is used for visualization
  void export_data_model_to_dir(const std::string& prefix, const bool only_for_viz = true) const;
  void export_data_model(const std::string& file_name, const bool only_for_viz) const;

  void to_data_model(Json::Value& root, const bool only_for_viz) const;

  // ---------------------- other ----------------------
  BNG::SpeciesContainer& get_all_species() { return bng_engine.get_all_species(); }
  const BNG::SpeciesContainer& get_all_species() const { return bng_engine.get_all_species(); }

  BNG::RxnContainer& get_all_rxns() { return bng_engine.get_all_rxns(); }
  const BNG::RxnContainer& get_all_rxns() const { return bng_engine.get_all_rxns(); }

  count_buffer_id_t create_count_buffer(const std::string file_name, const size_t buffer_size) {
    count_buffer_id_t id = count_buffers.size();
    count_buffers.push_back(CountBuffer(file_name, buffer_size));
    return id;
  }

  CountBuffer& get_count_buffer(const count_buffer_id_t id) {
    assert(id < count_buffers.size());
    return count_buffers[id];
  }

  const CountBuffer& get_count_buffer(const count_buffer_id_t id) const {
    assert(id < count_buffers.size());
    return count_buffers[id];
  }

  GeometryObject& get_geometry_object(const geometry_object_id_t id) {
    // TODO: there will be multiple places where geom object id and index are mixed
    return get_partition(0).get_geometry_object_by_id(id);
  }

  const GeometryObject& get_geometry_object(const geometry_object_id_t id) const {
    // TODO: there will be multiple places where geom object id and index are mixed
    return get_partition(0).get_geometry_object_by_id(id);
  }

  const Region& get_region(const region_id_t id) const {
    return get_partition(0).get_region_by_id(id);
  }

  const Region* find_region_by_name(const std::string& name) const {
    return get_partition(0).find_region_by_name(name);
  }

  static uint64_t determine_output_frequency(uint64_t iterations);

  bool check_for_overlapped_walls();

  void reset_unimol_rxn_times(const BNG::rxn_rule_id_t rxn_rule_id);

  // gives ownership of the event to this World object
  void add_unscheduled_count_event(MolOrRxnCountEvent* e) {
    unscheduled_count_events.push_back(e);
  }

  void flush_buffers();

private:
  // called in init_simulation
  void recompute_species_flags();

  void init_fpu();
  void init_counted_volumes();

  void initialization_to_data_model(Json::Value& mcell_node) const;

  void export_data_layout() const;

  std::string export_releases_to_bngl_seed_species(
      std::ostream& parameters, std::ostream& seed_species) const;

  std::string export_counts_to_bngl_observables(std::ostream& observables) const;
public:
  // single instance for the whole mcell simulator,
  // used as constants during simulation
  SimulationConfig config;

  BNG::BNGEngine bng_engine;

  SimulationStats stats;

  // owned by API::Model or references a global instance in case of MDL mode
  API::Callbacks& callbacks;

  SpeciesFlagsAnalyzer species_flags_analyzer;

  Scheduler scheduler;

  std::vector<MolOrRxnCountEvent*> unscheduled_count_events;

  uint64_t total_iterations; // number of iterations to simulate - move to Sim config

  rng_state rng; // single state for the random number generator

private:
  PartitionVector partitions;

  CountBufferVector count_buffers;

  // periodic check of used memory using timer
  MemoryLimitChecker memory_limit_checker;

  // global ID counters
  wall_id_t next_wall_id;
  region_id_t next_region_id;
  geometry_object_id_t next_geometry_object_id;

  // used by run_n_iterations to know whether the simulation was
  // already initialized
  bool simulation_initialized;

  // used to know whether we already reported final simulation stats
  bool simulation_ended;

  // buffers can be flushed only once
  bool buffers_flushed;

  // several variables to report simulation time
  timeval previous_progress_report_time;
  rusage sim_start_time;
  bool it1_start_time_set;
  rusage it1_start_time; // time when 1st iteration started

  // and to nicely report simulation progress
  uint64_t previous_iteration;
};

} // namespace mcell

#endif // SRC4_WORLD_H_

