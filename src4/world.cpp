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

#include <time.h>
#include <sys/time.h> // Linux include
#include <sys/resource.h> // Linux include
#include <fenv.h> // Linux include

//extern "C" {
#include "rng.h" // MCell 3
#include "logging.h"
//}

#include "world.h"
#include "end_simulation_event.h"
#include "defragmentation_event.h"

using namespace std;

const double USEC_IN_SEC = 1000000.0;

namespace MCell {

World::World()
  : current_iteration(0),
    iterations(0),
    seed_seq(0),
    next_wall_id(0),
    next_geometry_object_id(0)
{
  world_constants.partition_edge_length = PARTITION_EDGE_LENGTH_DEFAULT;
  world_constants.subpartitions_per_partition_dimension = SUBPARTITIONS_PER_PARTITION_DIMENSION_DEFAULT;
}


void World::init_fpu() {
#ifdef NDEBUG
  // we do not want to be making extra checks for division by zero
  // all places where such a case can occur is marked with comment POSSIBLE ZERO DIV
  fedisableexcept(FE_DIVBYZERO);

  static float_t a = 1;
  static float_t b = 0;
  if (a/b != INFINITY) {
    mcell_error("Error: division by zero is expected to return INFINITY but does not!");
    exit(1);
  }
#endif
}

void World::init_world_constants() {

  if (species.empty()) {
    // for developers: init_world_constants must be called after species and reactions conversion
    mcell_log("Error: there must be at lease one species!");
    exit(1);
  }

  // create map for fast reaction searches
  for (Reaction& r: reactions) {
    assert(r.reactants.size() == 1 || r.reactants.size() == 2); // only bimolecular reactions are supported now

    if (r.reactants.size() == 1) {
      // for now we only support only one outcome of a bimolecular reaction
      assert(unimolecular_reactions_map.count(r.reactants[0].species_id) == 0);
      unimolecular_reactions_map[r.reactants[0].species_id] = &r;
    }
    else {
      // check - for now we only support only one outcome of a bimolecular reaction
      if (bimolecular_reactions_map.count(r.reactants[0].species_id) != 0) {
        assert(bimolecular_reactions_map[r.reactants[0].species_id].count(r.reactants[1].species_id) == 0);
      }
      bimolecular_reactions_map[r.reactants[0].species_id][r.reactants[1].species_id] = &r;
      bimolecular_reactions_map[r.reactants[1].species_id][r.reactants[0].species_id] = &r;
    }
  }

  // just to make sure that we have an item for all the species
  for (Species& s: species) {
    bimolecular_reactions_map.insert( std::make_pair(s.species_id, SpeciesReactionMap()) );
  }
  assert(bimolecular_reactions_map.size() == species.size());

  world_constants.init(&unimolecular_reactions_map, &bimolecular_reactions_map, &species);
}


void World::init_simulation() {

  init_fpu();

  cout << "Partitions contain " <<  world_constants.subpartitions_per_partition_dimension << "^3 subvolumes.";
  assert(partitions.size() == 1 && "Initial partition must have been created, only 1 for now");
}


static uint64_t determine_output_frequency(uint64_t iterations) {
  uint64_t frequency;

  if (iterations < 10)
    frequency = 1;
  else if (iterations < 1000)
    frequency = 10;
  else if (iterations < 100000)
    frequency = 100;
  else if (iterations < 10000000)
    frequency = 1000;
  else if (iterations < 1000000000)
    frequency = 10000;
  else
    frequency = 100000;

  return frequency;
}


static double tousecs(timeval& t) {
  return (double)t.tv_sec * USEC_IN_SEC + (double)t.tv_usec;
}


static double tosecs(timeval& t) {
  return (double)t.tv_sec + (double)t.tv_usec/USEC_IN_SEC;
}


bool World::run_simulation(const bool dump_initial_state) {

  init_simulation(); // must be the first one

  if (dump_initial_state) {
    dump();
  }

  // create defragmentation events
  DefragmentationEvent* defragmentation_event = new DefragmentationEvent(this);
  defragmentation_event->event_time = DEFRAGMENTATION_PERIODICITY;
  defragmentation_event->periodicity_interval = DEFRAGMENTATION_PERIODICITY;
  scheduler.schedule_event(defragmentation_event);

  // create event that will terminate our simulation
  EndSimulationEvent* end_event = new EndSimulationEvent();
  end_event->event_time = iterations;
  scheduler.schedule_event(end_event);

  bool end_simulation = false;
  current_iteration = 0;
  uint64_t previous_iteration = 0;
  uint output_frequency = determine_output_frequency(iterations);
  timeval last_timing_time = {0, 0};

  cout << "Iterations: " << current_iteration << " of " << iterations << "\n";

  rusage sim_start_time;
  getrusage(RUSAGE_SELF, &sim_start_time);

  do {
#ifdef DEBUG_SCHEDULER
    cout << "Before it: " << current_iteration << ", time: " << time << "\n";
#endif

    // current_iteration corresponds to the number of executed time steps
    float_t time = scheduler.get_next_event_time();
    current_iteration = (uint64_t)time;

    // this is where events get executed
    scheduler.handle_next_event(end_simulation);

    // report progress
    if (current_iteration > previous_iteration) {

      if (current_iteration % output_frequency == 0) {
        cout << "Iterations: " << current_iteration << " of " << iterations;

        timeval curr_timing_time;
        gettimeofday(&curr_timing_time, NULL);
        if (last_timing_time.tv_usec > 0) {
          double time_diff = tousecs(curr_timing_time) - tousecs(last_timing_time);
          time_diff /= (double)output_frequency;
          cout << " (" << 1000000.0/time_diff << " iter/sec)";
        }
        last_timing_time = curr_timing_time;

        cout << "\n";
      }

      previous_iteration = current_iteration;
    }

#ifdef DEBUG_SCHEDULER
    cout << "After it: " << current_iteration << ", time: " << time << "\n";
#endif

  } while (!end_simulation);

  cout << "Iteration " << current_iteration << ", simulation finished successfully\n";

  simulation_stats.dump();

  // report final time
  rusage run_time;
  getrusage(RUSAGE_SELF, &run_time);
  cout << "Simulation CPU time = "
    << tosecs(run_time.ru_utime) - tosecs(sim_start_time.ru_utime) <<  "(user) and "
    << tosecs(run_time.ru_stime) - tosecs(sim_start_time.ru_stime) <<  "(system)\n";

  return true;
}


void World::dump() {
  world_constants.dump();
  // species
  Species::dump_array(species);

  // partitions
  for (Partition& p: partitions) {
    p.dump();
  }
}

} // namespace mcell

