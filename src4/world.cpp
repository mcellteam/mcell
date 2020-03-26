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

#include <fenv.h> // Linux include
#include <sys/resource.h> // Linux include

#include "rng.h" // MCell 3
#include "logging.h"


#include "world.h"
#include "defragmentation_event.h"

using namespace std;

const double USEC_IN_SEC = 1000000.0;

namespace MCell {

World::World()
  : iterations(0),
    seed_seq(0),
    next_wall_id(0),
    next_geometry_object_id(0),
    simulation_initialized(false),
    simulation_ended(false),
    previous_progress_report_time({0, 0}),
    previous_iteration(0),
    wall_hit_callback(nullptr),
    wall_hit_callback_clientdata(nullptr)
{
  config.partition_edge_length = PARTITION_EDGE_LENGTH_DEFAULT;
  config.subpartitions_per_partition_dimension = SUBPARTITIONS_PER_PARTITION_DIMENSION_DEFAULT;
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




void World::init_simulation() {
  assert(!simulation_initialized && "init_simulation must be called just once");

  if (bng_engine.all_species.get_count() == 0) {
    mcell_log("Error: there must be at lease one species!");
    exit(1);
  }

  // TODO: what do I need for initialization?
  config.init();
  stats.reset();

  init_fpu();

  cout << "Partitions contain " <<  config.subpartitions_per_partition_dimension << "^3 subvolumes.";
  assert(partitions.size() == 1 && "Initial partition must have been created, only 1 is allowed for now");

  // create defragmentation events
  DefragmentationEvent* defragmentation_event = new DefragmentationEvent(this);
  defragmentation_event->event_time = DEFRAGMENTATION_PERIODICITY;
  defragmentation_event->periodicity_interval = DEFRAGMENTATION_PERIODICITY;
  scheduler.schedule_event(defragmentation_event);

  // initialize timing
  previous_progress_report_time = {0, 0};

  rusage sim_start_time;
  getrusage(RUSAGE_SELF, &sim_start_time);

  // iteration counter to report progress
  previous_iteration = 0;

  simulation_initialized = true;
}


void World::run_n_iterations(const uint64_t num_iterations, const uint64_t output_frequency, const bool terminate_last_iteration_after_viz_output) {

  if (!simulation_initialized) {
    init_simulation();
  }

  uint64_t& current_iteration = stats.get_current_iteration();

  if (current_iteration == 0) {
    cout << "Iterations: " << current_iteration << " of " << iterations << "\n";
  }

  uint64_t this_run_first_iteration = current_iteration;

  do {
    // current_iteration corresponds to the number of executed time steps
    float_t time = scheduler.get_next_event_time();
    current_iteration = (uint64_t)time;

#ifdef DEBUG_SCHEDULER
    cout << "Before it: " << current_iteration << ", time: " << time << "\n";
#endif

    // terminate simulation if we executed the right number of iterations
    if (current_iteration >= this_run_first_iteration + num_iterations) {
      break;
    }

    // this is where events get executed
    EventExecutionInfo event_info = scheduler.handle_next_event();

    // report progress
    if (current_iteration > previous_iteration) {

      if (current_iteration % output_frequency == 0) {
        cout << "Iterations: " << current_iteration << " of " << iterations;

        timeval current_progress_report_time;
        gettimeofday(&current_progress_report_time, NULL);
        if (previous_progress_report_time.tv_usec > 0) {
          double time_diff = tousecs(current_progress_report_time) - tousecs(previous_progress_report_time);
          time_diff /= (double)output_frequency;
          cout << " (" << 1000000.0/time_diff << " iter/sec)";
        }
        previous_progress_report_time = current_progress_report_time;

        cout << "\n";
      }

      previous_iteration = current_iteration;
    }

#ifdef DEBUG_SCHEDULER
    cout << "After it: " << current_iteration << ", time: " << time << "\n";
#endif

    // also terminate if this was the last iteration and the event was viz output
    if (terminate_last_iteration_after_viz_output &&
       current_iteration == this_run_first_iteration + num_iterations - 1 &&
       event_info.type_index >= EVENT_TYPE_INDEX_VIZ_OUTPUT
    ) {
      break;
    }

  } while (true); // terminated when the nr. of iterations is reached

#ifndef NDEBUG
  // flush everything, we want the output to be mixed with Python in the right ordering
  cout.flush();
  cerr.flush();
  fflush(stdout);
  fflush(stderr);
#endif
}


void World::end_simulation() {
  if (simulation_ended) {
    // already called, do nothing
    return;
  }

  cout << "Iteration " << stats.get_current_iteration() << ", simulation finished successfully\n";

  stats.dump();

  // report final time
  rusage run_time;
  getrusage(RUSAGE_SELF, &run_time);
  cout << "Simulation CPU time = "
    << tosecs(run_time.ru_utime) - tosecs(sim_start_time.ru_utime) <<  "(user) and "
    << tosecs(run_time.ru_stime) - tosecs(sim_start_time.ru_stime) <<  "(system)\n";

  simulation_ended = true;
}


void World::run_simulation(const bool dump_initial_state) {

  // do initialization, also insert
  // defragmentation and end simulation event
  init_simulation();

  if (dump_initial_state) {
    dump();
  }

  uint output_frequency = determine_output_frequency(iterations);

  // simulating 1000 iterations means to simulate iterations 0 .. 1000
  run_n_iterations(iterations + 1, output_frequency, true);

  end_simulation();
}


void World::dump() {
  stats.dump();

  bng_engine.all_species.dump();
  all_reactions.dump();

  // partitions
  for (Partition& p: partitions) {
    p.dump();
  }
}

} // namespace mcell

