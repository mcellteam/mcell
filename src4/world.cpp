/******************************************************************************
 *
 * Copyright (C) 2019-2021 by
 * The Salk Institute for Biological Studies
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

#ifndef _MSC_VER
#include <sys/time.h>
#include <sys/resource.h>
#include <sys/resource.h> // Linux include
#endif
#include <signal.h>

#include <fstream>

#include "rng.h" // MCell 3
#include "logging.h"
#include "util.h"

#include "world.h"
#include "viz_output_event.h"
#include "defragmentation_event.h"
#include "rxn_class_cleanup_event.h"
#include "species_cleanup_event.h"
#include "sort_mols_by_subpart_event.h"
#include "release_event.h"
#include "mol_or_rxn_count_event.h"
#include "datamodel_defines.h"
#include "bng_data_to_datamodel_converter.h"
#include "diffuse_react_event.h"
#include "run_n_iterations_end_event.h"
#include "custom_function_call_event.h"
#include "mol_order_shuffle_event.h"
#include "vtk_utils.h"

#include "api/mol_wall_hit_info.h"
#include "api/geometry_object.h"
#include "api/model.h"
#include "generated/gen_constants.h"

#include "bng/filesystem_utils.h"

using namespace std;

const double USEC_IN_SEC = 1000000.0;

namespace MCell {

static double tousecs(timeval& t) {
  return (double)t.tv_sec * USEC_IN_SEC + (double)t.tv_usec;
}


static double tosecs(timeval& t) {
  return (double)t.tv_sec + (double)t.tv_usec/USEC_IN_SEC;
}


static void print_periodic_stats_func(double time, World* world) {
  release_assert(world != nullptr);
  world->print_periodic_stats();
}


World::World(API::Callbacks& callbacks_)
  : bng_engine(config),
    callbacks(callbacks_),
    total_iterations(0),
    next_wall_id(0),
    next_region_id(0),
    next_geometry_object_id(0),
    simulation_initialized(false),
    run_n_iterations_terminated_with_checkpoint(false),
    simulation_ended(false),
    buffers_flushed(false),
    it1_start_time_set(false),
    previous_iteration(0),
    signaled_checkpoint_signo(API::SIGNO_NOT_SIGNALED),
    signaled_checkpoint_model(nullptr)
{
  config.partition_edge_length = FLT_INVALID;
  config.num_subparts_per_partition_edge = SUBPARTITIONS_PER_PARTITION_DIMENSION_DEFAULT;

  // although the same thing is called in init_simulation, not reseting it causes weird valdrind reports on
  // uninitialized variable
  reset_rusage(&sim_start_time);
}


World::~World() {
  if (!buffers_flushed) {
    flush_and_close_buffers();
  }

  for (MolOrRxnCountEvent* e: unscheduled_count_events) {
    delete e;
  }
}


void World::init_fpu() {
  // default is FE_TONEAREST but let's set it to be sure
  fesetround(FE_TONEAREST);
}


uint64_t World::determine_output_frequency(uint64_t iterations) {
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


void World::create_initial_surface_region_release_event() {
  ReleaseEvent* rel_event = new ReleaseEvent(this);
  rel_event->event_time = 0;
  rel_event->release_site_name = "initial surface releases";
  rel_event->release_shape = ReleaseShape::INITIAL_SURF_REGION;
  rel_event->update_event_time_for_next_scheduled_time();
  scheduler.schedule_event(rel_event);
}


void World::init_counted_volumes() {
  assert(partitions.size() == 1);

  bool ok = VtkUtils::initialize_counted_volumes(this, config.has_intersecting_counted_objects);
  if (!ok) {
    mcell_error("Processing of counted volumes failed, terminating.");
  }

  partitions[PARTITION_ID_INITIAL].initialize_all_waypoints();
}


static double get_event_start_time(const double start_time, const double periodicity) {
  if (periodicity == 0) {
    return 0;
  }
  else {
    return floor_to_multiple_f(start_time + periodicity, periodicity);
  }
}


void World::schedule_checkpoint_event(
    const uint64_t iteration, const bool continue_simulation, const API::CheckpointSaveEventContext& ctx) {

  CustomFunctionCallEvent<API::CheckpointSaveEventContext>* checkpoint_event =
      new CustomFunctionCallEvent<API::CheckpointSaveEventContext>(
          API::save_checkpoint_func, ctx, EVENT_TYPE_INDEX_CALL_START_ITERATION_CHECKPOINT);

  if (iteration == 0) {
    // schedule for the closest iteration while correctly maintaining order of events in the queue
    checkpoint_event->event_time = TIME_INVALID;
  }
  else {
    checkpoint_event->event_time = iteration;
  }
  checkpoint_event->periodicity_interval = 0; // only once
  checkpoint_event->return_from_run_n_iterations = !continue_simulation;

  // safely schedule
  // not really needed to do safely when called due to a signal handler because only a flag
  // is set and the handling is done synchronously, but required when schedule_checkpoint is called e.g. from
  // a timer
  scheduler.schedule_event_asynchronously(checkpoint_event);
}


void World::check_checkpointing_signal() {
  if (signaled_checkpoint_signo == API::SIGNO_NOT_SIGNALED) {
    return;
  }

  bool continue_simulation = false;

  // printout for each model instance
#ifndef _WIN32
  if (signaled_checkpoint_signo == SIGUSR1) {
    cout << "User signal SIGUSR1 detected, scheduling a checkpoint and continuing simulation.\n";
    continue_simulation = true;
  }
  else if (signaled_checkpoint_signo == SIGUSR2) {
    cout << "User signal SIGUSR2 detected, scheduling a checkpoint and terminating simulation afterwards.\n";
    continue_simulation = false;
  }
  else
#endif
  if (signaled_checkpoint_signo == SIGALRM) {
    cout << "Signal SIGALRM detected - periodic or time limit elapsed, scheduling a checkpoint ";
    if (config.continue_after_sigalrm) {
      cout << "and continuing simulation.\n";
      continue_simulation = true;
    }
    else {
      cout << "and terminating simulation afterwards.\n";
      continue_simulation = false;
    }
  }
  else {
    cout << "Unexpected signal " << signaled_checkpoint_signo << " received, fatal error.\n";
    release_assert(false);
  }

  API::CheckpointSaveEventContext ctx;
  release_assert(signaled_checkpoint_model != nullptr);
  ctx.model = signaled_checkpoint_model;
  ctx.dir_prefix = config.get_default_checkpoint_dir_prefix();
  ctx.append_it_to_dir = true;

  schedule_checkpoint_event(0, continue_simulation, ctx);

  // reset flag and pointer to model because we don't need it anymore
  signaled_checkpoint_signo = API::SIGNO_NOT_SIGNALED;
  signaled_checkpoint_model = nullptr;
}


void World::init_simulation(const double start_time) {

  release_assert((int)start_time == start_time && "Iterations start time must be an integer.");

  // TODO: check these messages in testsuite
#ifdef MCELL3_4_ALWAYS_SORT_MOLS_BY_TIME_AND_ID
  mcell_log("!!! WARNING: Event sorting according to time and id was enabled for debugging, testing won't pass.");
#endif
#ifdef MCELL4_DO_NOT_REUSE_REACTANT
  mcell_log("!!! WARNING: Reactant reuse is disabled, testing won't pass.");
#endif
#ifdef MCELL4_SORT_RXN_PRODUCTS_BY_NAME
  mcell_log("!!! WARNING: Standard product sorting is disabled, testing won't pass.");
#endif

  if (DUMP4_PRECISION != DUMP4_PRECISION_DEFAULT) {
    cout.precision(DUMP4_PRECISION);
  }

  assert(!simulation_initialized && "init_simulation must be called just once");

  if (get_all_species().get_count() == 0) {
    mcell_log("Error: there must be at least one species!");
    exit(1);
  }

  stats.reset(false);

  init_fpu();

  init_counted_volumes();

  cout <<
      "Partition contains " <<  config.num_subparts_per_partition_edge << "^3 subpartitions, " <<
      "subpartition size is " << config.subpart_edge_length * config.length_unit << " microns.\n";
  assert(partitions.size() == 1 && "Initial partition must have been created, only 1 is allowed for now");

  // create event that diffuses molecules
  DiffuseReactEvent* event = new DiffuseReactEvent(this);
  event->event_time = start_time;
  scheduler.schedule_event(event);

  // create defragmentation events
  DefragmentationEvent* defragmentation_event = new DefragmentationEvent(this);
  defragmentation_event->event_time = get_event_start_time(start_time, DEFRAGMENTATION_PERIODICITY);
  defragmentation_event->periodicity_interval = DEFRAGMENTATION_PERIODICITY;
  scheduler.schedule_event(defragmentation_event);

  // create rxn class cleanup events
  if (config.rxn_class_cleanup_periodicity >= 1) {
    RxnClassCleanupEvent* rxn_class_cleanup_event = new RxnClassCleanupEvent(this);
    rxn_class_cleanup_event->event_time = get_event_start_time(start_time, config.rxn_class_cleanup_periodicity);
    rxn_class_cleanup_event->periodicity_interval = config.rxn_class_cleanup_periodicity;
    scheduler.schedule_event(rxn_class_cleanup_event);
  }

  if (config.species_cleanup_periodicity >= 1) {
    SpeciesCleanupEvent* species_cleanup_event = new SpeciesCleanupEvent(this);
    species_cleanup_event->event_time = get_event_start_time(start_time, config.species_cleanup_periodicity);
    species_cleanup_event->periodicity_interval = config.species_cleanup_periodicity;
    scheduler.schedule_event(species_cleanup_event);
  }

  if (config.molecules_order_random_shuffle_periodicity >= 1) {
    MolOrderShuffleEvent* mol_order_shuffle_event = new MolOrderShuffleEvent(this);
    mol_order_shuffle_event->event_time = get_event_start_time(start_time, config.molecules_order_random_shuffle_periodicity);
    mol_order_shuffle_event->periodicity_interval = config.molecules_order_random_shuffle_periodicity;
    scheduler.schedule_event(mol_order_shuffle_event);
  }

  // create subpart sorting events
  if (config.sort_mols_by_subpart) {
    SortMolsBySubpartEvent* sort_event = new SortMolsBySubpartEvent(this);
    if (start_time == 0) {
      sort_event->event_time = TIME_SIMULATION_START;
    }
    else {
      sort_event->event_time = get_event_start_time(start_time, SORT_MOLS_BY_SUBPART_PERIODICITY);
    }
    sort_event->periodicity_interval = SORT_MOLS_BY_SUBPART_PERIODICITY;
    scheduler.schedule_event(sort_event);
  }

  // simulation statistics, mostly for development purposes
  if (config.simulation_stats_every_n_iterations > 0) {
    CustomFunctionCallEvent<World*>* stats_event =
        new CustomFunctionCallEvent<World*>(print_periodic_stats_func, this);
    stats_event->event_time = start_time;
    stats_event->periodicity_interval = config.simulation_stats_every_n_iterations;
    scheduler.schedule_event(stats_event);
  }

  // initialize timing
  previous_progress_report_time = chrono::steady_clock::now();
  previous_buffer_flush_time = previous_progress_report_time;
  reset_rusage(&sim_start_time);
  getrusage(RUSAGE_SELF, &sim_start_time);

  // iteration counter to report progress
  previous_iteration = 0;

  if (!config.use_expanded_list) {
    cout <<
        "Warning: configuration 'use_expanded_list' (ACCURATE_3D_REACTIONS) set to false is compatible "
        "with MCell3 only in simple cases without volume reactions and usually MCell4 produces different results. "
        "Search for potential reactions in MCell3 is always based on reaction radius, "
        "the solution in MCell3 can be highly imprecise because it considers subvolume/subpartition boundaries "
        "when computing reaction probability factor in exact-disk.\n";
  }

  // start memory check timer
  memory_limit_checker.start_timed_check(this, config.memory_limit_gb);

  if (stats.get_current_iteration() == 0) {
    // not starting from a checkpoint
    config.initialize_run_report_file();
    BNG::append_to_report(config.get_run_report_file_name(), "Simulation started ");
  }
  else {
    cout << "Iterations: " << stats.get_current_iteration() << " of " << total_iterations << " (resuming a checkpoint)\n";
    BNG::append_to_report(config.get_run_report_file_name(), "Simulation resumed from a checkpoint ");
  }
  BNG::append_to_report(config.get_run_report_file_name(),
      "at iteration " + to_string(stats.get_current_iteration()) +
      " and time " + BNG::get_current_date_time() + ".\n");

  simulation_initialized = true;
}


uint64_t World::time_to_iteration(const double time) {
  return (uint64_t)(time - config.get_simulation_start_time()) + config.initial_iteration;
}


uint64_t World::run_n_iterations(const uint64_t num_iterations, const bool terminate_last_iteration_after_viz_output) {

  release_assert(simulation_initialized);

  run_n_iterations_terminated_with_checkpoint = false;
  uint64_t output_frequency = determine_output_frequency(total_iterations);
  uint64_t this_run_first_iteration = stats.get_current_iteration();
  uint64_t& current_iteration = stats.get_current_iteration();

  if (current_iteration == 0 && config.iteration_report) {
    cout << "Iterations: " << current_iteration << " of " << total_iterations << "\n";
  }

  // information for scheduling of RunNIterationsEndEvent
  bool run_n_iters_end_event_created = false;
  uint64_t iteration_to_create_run_n_iters_end_event;
  if (num_iterations <= ITERATIONS_BEFORE_RUN_N_ITERATIONS_END_EVENT) {
    // schedule right away
    iteration_to_create_run_n_iters_end_event = this_run_first_iteration;
  }
  else {
    // schedule later
    iteration_to_create_run_n_iters_end_event =
        this_run_first_iteration + num_iterations - ITERATIONS_BEFORE_RUN_N_ITERATIONS_END_EVENT;
  }

  do {
    check_checkpointing_signal();

    if (!run_n_iters_end_event_created && iteration_to_create_run_n_iters_end_event == current_iteration) {
      // create event that is used to check whether simulation should end right after the last viz output,
      // also serves as a simulation barrier to not to do diffusion after this point in time
      // must not be created right away when run_n_iterations is called because this might require too much memory
      RunNIterationsEndEvent* run_n_iterations_end_event = new RunNIterationsEndEvent();
      run_n_iterations_end_event->event_time = this_run_first_iteration + num_iterations;
      run_n_iterations_end_event->periodicity_interval = 0; // these markers are inserted into every time step
      scheduler.schedule_event(run_n_iterations_end_event);
      run_n_iters_end_event_created = true;
    }

    // current_iteration corresponds to the number of executed time steps
    double time = scheduler.get_next_event_time();

    // convert time to iteration
    current_iteration = time_to_iteration(time);

    if (current_iteration == 1 && previous_iteration == 0) {
      it1_start_time_set = true;
      reset_rusage(&it1_start_time);
      getrusage(RUSAGE_SELF, &it1_start_time);
    }

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
        auto current_time = std::chrono::steady_clock::now();

        if (config.iteration_report) {
          cout << "Iterations: " << current_iteration << " of " << total_iterations;

          double iters_per_sec =
            1000000.0 /
            ((chrono::duration_cast<chrono::microseconds>(current_time - previous_progress_report_time).count() /
             (double)output_frequency));
          cout << " (" << iters_per_sec << " iter/sec)";
          previous_progress_report_time = current_time;

          cout << " " << bng_engine.get_stats_report();

          cout << "\n";
          cout.flush(); // flush is required so that CellBlender can display progress
        }

        // should we also flush count buffers?
        double time_diff = chrono::duration_cast<chrono::seconds>(current_time - previous_buffer_flush_time).count();
        if (time_diff >= COUNT_BUFFER_FLUSH_SECONDS) {
          flush_buffers();
          previous_buffer_flush_time = current_time;
        }
      }

      previous_iteration = current_iteration;
    }

#ifdef DEBUG_SCHEDULER
    cout << "After it: " << current_iteration << ", time: " << time << "\n";
#endif

    // also terminate if this was the last iteration and we hit an event that represents a check for the
    // end of the simulation
    if (
        (terminate_last_iteration_after_viz_output &&
         event_info.type_index == EVENT_TYPE_INDEX_SIMULATION_END_CHECK
        ) ||
        event_info.return_from_run_iterations
    ) {
      assert(
          event_info.return_from_run_iterations ||
          current_iteration == this_run_first_iteration + num_iterations - 1);

      if (event_info.type_index == EVENT_TYPE_INDEX_CALL_START_ITERATION_CHECKPOINT) {
        run_n_iterations_terminated_with_checkpoint = true;
      }

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

  return current_iteration - this_run_first_iteration;
}


count_buffer_id_t World::create_dat_count_buffer(
    const std::string file_name, const size_t buffer_size, const bool open_for_append) {
  count_buffer_id_t id = count_buffers.size();
  std::vector<std::string> column_names = { file_name }; // name is not used when .dat format is used
  count_buffers.push_back(
      CountBuffer(CountOutputFormat::DAT, file_name, column_names, buffer_size, open_for_append));
  count_buffers.back().open();
  return id;
}


count_buffer_id_t World::create_gdat_count_buffer(
    const std::string file_name, const std::vector<std::string>& column_names,
    const size_t buffer_size, const bool open_for_append) {
  count_buffer_id_t id = count_buffers.size();
  count_buffers.push_back(
      CountBuffer(CountOutputFormat::GDAT, file_name, column_names, buffer_size, open_for_append));
  count_buffers.back().open();
  return id;
}


void World::flush_buffers() {
  // only flush count buffers
  for (CountBuffer& b: count_buffers) {
    b.flush();
  }
}


void World::flush_and_close_buffers() {
  assert(!buffers_flushed && "Buffers can be flushed only once");

  // flush and close count buffers
  for (CountBuffer& b: count_buffers) {
    b.flush_and_close();
  }
  buffers_flushed = true;
}


// prints message (appends newline), flushes buffers, and terminates
void World::fatal_error(const std::string& msg) {
  errs() << msg << "\n";
  flush_and_close_buffers();
  exit(1);
}


void World::end_simulation(const bool print_final_report) {
  // we do not want to check memory anymore
  memory_limit_checker.stop_timed_check();

  if (simulation_ended) {
    // already called, do nothing
    return;
  }

  // - execute all events up to the last viz output
  //   to produce viz output and counts for the last iteration
  // - must not be done if we ended with checkpoint
  if (!run_n_iterations_terminated_with_checkpoint) {
    run_n_iterations(1, true);
  }

  flush_and_close_buffers();

  if (print_final_report) {
    cout << "Iteration " << stats.get_current_iteration() << ", simulation finished successfully";
    if (run_n_iterations_terminated_with_checkpoint) {
      cout << ", terminated with a checkpoint";
    }
    cout << "\n";

    stats.print_report();

    // report final time
    rusage run_time;
    reset_rusage(&run_time);
    getrusage(RUSAGE_SELF, &run_time);
    cout << "Simulation CPU time = "
      << tosecs(run_time.ru_utime) - tosecs(sim_start_time.ru_utime) <<  " (user) and "
      << tosecs(run_time.ru_stime) - tosecs(sim_start_time.ru_stime) <<  " (system)\n";
    cout << "Simulation CPU time without iteration 0 = "
      << tosecs(run_time.ru_utime) - tosecs(it1_start_time.ru_utime) <<  " (user) and "
      << tosecs(run_time.ru_stime) - tosecs(it1_start_time.ru_stime) <<  " (system)\n";

    // and warnings
    bng_engine.get_config().print_final_warnings();
  }

  BNG::append_to_report(config.get_run_report_file_name(),
      "Simulation ended at iteration " + to_string(stats.get_current_iteration()) +
      " and time " + BNG::get_current_date_time() + ", " +
      ((run_n_iterations_terminated_with_checkpoint) ?
          "terminated due to checkpoint.\n" :
          "all iterations were finished.\nFINISHED\n"));

  simulation_ended = true;
}


void World::init_and_run_simulation(const bool dump_initial_state, const bool dump_with_geometry) {

  // do initialization, also insert
  // defragmentation and end simulation event
  init_simulation(TIME_SIMULATION_START);

  if (dump_initial_state) {
    dump(dump_with_geometry);
  }

  run_n_iterations(total_iterations, true);

  // runs one more iteration but only up to the last viz output
  end_simulation(true);
}


void World::print_periodic_stats() const {
  cout << "--- Periodic Stats ---\n";
  cout << "World: stats.current_iteration = " << stats.get_current_iteration() << "\n";
  scheduler.print_periodic_stats();
  bng_engine.print_periodic_stats();
  get_partition(PARTITION_ID_INITIAL).print_periodic_stats();
  cout << "Memory: " << get_mem_usage() << " kB\n";
  cout << "----------------------\n";
}


void World::dump(const bool with_geometry) {
  config.dump();
  stats.print_report();

  bng_engine.get_data().dump();
  get_all_species().dump();
  get_all_rxns().dump(true);

  // partitions
  for (Partition& p: partitions) {
    p.dump(with_geometry);
  }

  scheduler.dump();
}


bool World::check_for_overlapped_walls() {
  /* pick up a random vector */
  Vec3 rand_vec;
  rand_vec.x = rng_dbl(&rng);
  rand_vec.y = rng_dbl(&rng);
  rand_vec.z = rng_dbl(&rng);

  for (const Partition& p: partitions) {
    bool ok = p.check_for_overlapped_walls(rand_vec);
    if (!ok) {
      return false;
    }
  }
  return true;
}


void World::reset_unimol_rxn_times(const BNG::rxn_rule_id_t rxn_rule_id) {
  // get all affected species
  const BNG::RxnRule* rxn = get_all_rxns().get(rxn_rule_id);
  assert(rxn->is_unimol());

  // and then reset unimol time for each molecule of that species
  for (Partition& p: partitions) {
    for (Molecule& m: p.get_molecules()) {

      assert((rxn->species_applicable_as_any_reactant.count(m.species_id) != 0 ||
          rxn->species_not_applicable_as_any_reactant.count(m.species_id) != 0) &&
          "This unimol rxn must have been analyzed for this species because the species are instantiated");

      if (rxn->species_applicable_as_any_reactant.count(m.species_id) != 0) {
        // new unimol time will be computed when the molecule is diffused
        // the next time (we cannot change it right away for molecules that
        // have longer timestep anyway because they were already diffused to the future)

        // is this a nondiffusible molecule? - these don't have to be not scheduled and therefore
        // their unimol rxn time does not need to be recomputed
        const BNG::Species& s = get_all_species().get(m.species_id);
        if (!s.can_diffuse() && (m.diffusion_time == TIME_FOREVER || m.diffusion_time == m.unimol_rx_time)) {
          // we cannot change the unimol rate if we are currently diffusing, but
          // let's update it as soon as possible
          m.diffusion_time = scheduler.get_next_event_time();
        }

        m.unimol_rx_time = TIME_INVALID;
        m.clear_flag(MOLECULE_FLAG_RESCHEDULE_UNIMOL_RXN_ON_NEXT_RXN_RATE_UPDATE);
        m.set_flag(MOLECULE_FLAG_SCHEDULE_UNIMOL_RXN);
      }
    }
  }
}


void World::export_geometry_to_obj(const std::string& files_prefix) const {

  // create directories if needed
  FSUtils::make_dir_for_file_w_multiple_attempts(files_prefix);

  // only one partition now
  const Partition& p = get_partition(PARTITION_ID_INITIAL);

  VtkUtils::export_geometry_objects_to_obj(this, p.get_geometry_objects(), files_prefix);
}


void World::export_data_model_to_dir(const std::string& prefix, const bool only_for_viz) const {
  // prefix should be the same directory that is used for viz_output,
  // e.g. ./viz_data/seed_0001/Scene

  stringstream path;
  path <<
      prefix << ".data_model." <<
      VizOutputEvent::iterations_to_string(stats.get_current_iteration(), total_iterations) <<
      ".json";

  // create directories if needed
  FSUtils::make_dir_for_file_w_multiple_attempts(path.str());

  export_data_model(path.str(), only_for_viz);
}


void World::export_data_model(const std::string& file_name, const bool only_for_viz) const {

  Json::Value root;
  to_data_model(root, only_for_viz);

  Json::StreamWriterBuilder wbuilder;
  wbuilder["indentation"] = " ";
  wbuilder.settings_["precision"] = 15; // this is the precision that is used by mdl_to_data_model.py script
  wbuilder.settings_["precisionType"] = "significant";
  std::string document = Json::writeString(wbuilder, root);

  // write result into a file
  ofstream res_file(file_name);
  if (res_file.is_open())
  {
    res_file << document;
    res_file.close();
  }
  else {
    cout << "Unable to open file " << file_name << " for writing.\n";
  }

  // also, if rxn_output was enabled and in the first iteration, generate data_layout.json file
  // in the current directory
  if (stats.get_current_iteration() == 0) {
    export_data_layout();
  }
}


void World::to_data_model(Json::Value& root, const bool only_for_viz) const {
  Json::Value& mcell = root[KEY_MCELL];

  mcell[KEY_CELLBLENDER_VERSION] = VALUE_CELLBLENDER_VERSION;
  DMUtils::add_version(mcell, VER_DM_2017_06_23_1300);

  initialization_to_data_model(mcell);

  // generate geometry information

  // first create empty model_objects section (may be filled-in by
  // GeometryObject::to_data_model_as_model_object
  Json::Value& model_objects = mcell[KEY_MODEL_OBJECTS];
  DMUtils::add_version(model_objects, VER_DM_2018_01_11_1330);
  Json::Value& model_object_list = model_objects[KEY_MODEL_OBJECT_LIST];
  model_object_list = Json::Value(Json::arrayValue);

  // then dump all partition data
  set<rgba_t> used_colors;
  bool first = true;
  for (const Partition& p: partitions) {
    p.to_data_model(mcell, used_colors);
  }

  if (!only_for_viz) {
    // base information for reaction_data_output must be set even when there are no such events
    Json::Value& reaction_data_output = mcell[KEY_REACTION_DATA_OUTPUT];
    DMUtils::add_version(reaction_data_output, VER_DM_2016_03_15_1800);
    reaction_data_output[KEY_PLOT_LAYOUT] = " ";
    reaction_data_output[KEY_PLOT_LEGEND] = "0";
    reaction_data_output[KEY_MOL_COLORS] = false;
    reaction_data_output[KEY_ALWAYS_GENERATE] = true;
    reaction_data_output[KEY_OUTPUT_BUF_SIZE] = "";
    reaction_data_output[KEY_RXN_STEP] = "";
    reaction_data_output[KEY_COMBINE_SEEDS] = true;
    Json::Value& reaction_output_list = reaction_data_output[KEY_REACTION_OUTPUT_LIST];
    reaction_output_list = Json::Value(Json::arrayValue); // empty array

    scheduler.to_data_model(mcell);

  }

  // generate species and rxn info
  BngDataToDatamodelConverter bng_converter;
  bng_converter.to_data_model(this, mcell, only_for_viz);

  // add other default values, might need to generate this better
  Json::Value& materials = mcell[KEY_MATERIALS];
  Json::Value& material_dict = materials[KEY_MATERIAL_DICT];

  if (used_colors.empty()) {
    // need at least one material (unused)
    Json::Value& membrane = material_dict[KEY_VALUE_MEMBRANE];
    Json::Value& diffuse_color = membrane[KEY_DIFFUSE_COLOR];
    diffuse_color[KEY_R] = DEFAULT_OBJECT_COLOR_COMPONENT;
    diffuse_color[KEY_G] = DEFAULT_OBJECT_COLOR_COMPONENT;
    diffuse_color[KEY_B] = DEFAULT_OBJECT_COLOR_COMPONENT;
    diffuse_color[KEY_A] = DEFAULT_OBJECT_ALPHA;
  }
  else {
    for (rgba_t color: used_colors) {
      string name = DMUtils::color_to_mat_name(color);
      Json::Value& mat = material_dict[name];
      Json::Value& diffuse_color = mat[KEY_DIFFUSE_COLOR];

      double r, g, b, a;
      Geometry::rgba_to_components(color, r, g, b, a);
      diffuse_color[KEY_R] = r;
      diffuse_color[KEY_G] = g;
      diffuse_color[KEY_B] = b;
      diffuse_color[KEY_A] = a;
    }
  }

  // diverse settings not read from the mcell4 state

  // --- mol_viz ---
  Json::Value& mol_viz = mcell[KEY_MOL_VIZ];
  DMUtils::add_version(mol_viz, VER_DM_2015_04_13_1700);
  mol_viz[KEY_MANUAL_SELECT_VIZ_DIR] = false;
  mol_viz[KEY_FILE_START_INDEX] = 0;
  mol_viz[KEY_SEED_LIST] = Json::Value(Json::arrayValue); // empty array

  Json::Value& color_list = mol_viz[KEY_COLOR_LIST];
  DMUtils::append_triplet(color_list, 0.8, 0.0, 0.0);
  DMUtils::append_triplet(color_list, 0.0, 0.8, 0.0);
  DMUtils::append_triplet(color_list, 0.0, 0.0, 0.8);
  DMUtils::append_triplet(color_list, 0.0, 0.8, 0.8);
  DMUtils::append_triplet(color_list, 0.8, 0.0, 0.8);
  DMUtils::append_triplet(color_list, 0.8, 0.8, 0.0);
  DMUtils::append_triplet(color_list, 1.0, 1.0, 1.0);
  DMUtils::append_triplet(color_list, 0.0, 0.0, 0.0);

  mol_viz[KEY_ACTIVE_SEED_INDEX] = 0;
  mol_viz[KEY_FILE_INDEX] = 959; // don't know what this means
  mol_viz[KEY_FILE_NUM] = 1001; // don't know what this means
  mol_viz[KEY_VIZ_ENABLE] = true;
  mol_viz[KEY_FILE_NAME] = "";
  mol_viz[KEY_COLOR_INDEX] = 0;
  mol_viz[KEY_RENDER_AND_SAVE] = false;
  mol_viz[KEY_FILE_STEP_INDEX] = 1; // don't know what this means
  mol_viz[KEY_FILE_STOP_INDEX] = 1000; // don't know what this means
  mol_viz[KEY_FILE_DIR] = ""; // does this need to be set?

  mol_viz[KEY_VIZ_LIST] = Json::Value(Json::arrayValue); // empty array

  // --- parameter_system ---
  Json::Value& parameter_system = mcell[KEY_PARAMETER_SYSTEM];
  parameter_system[KEY_MODEL_PARAMETERS] = Json::Value(Json::ValueType::objectValue); // empty dict

  // --- scripting ---
  Json::Value& scripting = mcell[KEY_SCRIPTING];
  DMUtils::add_version(scripting, VER_DM_2017_11_30_1830);
  scripting[KEY_SCRIPTING_LIST] = Json::Value(Json::arrayValue);
  scripting[KEY_SCRIPT_TEXTS] = Json::Value(Json::objectValue);
  scripting[KEY_DM_INTERNAL_FILE_NAME] = "";
  scripting[KEY_FORCE_PROPERTY_UPDATE] = true;
  scripting[KEY_DM_EXTERNAL_FILE_NAME] = "";
  scripting[KEY_IGNORE_CELLBLENDER_DATA] = false;

  // --- simulation_control ---
  Json::Value& simulation_control = mcell[KEY_SIMULATION_CONTROL];
  simulation_control[KEY_EXPORT_FORMAT] = VALUE_MCELL_MDL_MODULAR;

  mcell[KEY_MODEL_LANGUAGE] = VALUE_MCELL4;

  Json::Value& blender_version = mcell[KEY_BLENDER_VERSION];
  blender_version.append(Json::Value(BLENDER_VERSION[0]));
  blender_version.append(Json::Value(BLENDER_VERSION[1]));
  blender_version.append(Json::Value(BLENDER_VERSION[2]));
}


void World::initialization_to_data_model(Json::Value& mcell_node) const {
  // only setting defaults for now, most of these values are not used in mcell4

  // --- initialization ---
  Json::Value& initialization = mcell_node[KEY_INITIALIZATION];
  DMUtils::add_version(initialization, VER_DM_2017_11_18_0130);

  // time step will most probably use rounded values, therefore we don't have to use full precision here
  initialization[KEY_TIME_STEP] = f_to_str(config.time_unit, 8);
  initialization[KEY_ITERATIONS] = to_string(total_iterations);

  if (!cmp_eq(config.rx_radius_3d, config.get_default_rx_radius_3d())) {
    initialization[KEY_INTERACTION_RADIUS] = f_to_str(config.rx_radius_3d * config.length_unit);
  }
  else {
    // keep default value
    initialization[KEY_INTERACTION_RADIUS] = "";
  }

  initialization[KEY_ACCURATE_3D_REACTIONS] = true;
  initialization[KEY_RADIAL_SUBDIVISIONS] = "";
  initialization[KEY_RADIAL_DIRECTIONS] = "";
  initialization[KEY_CENTER_MOLECULES_ON_GRID] = !config.randomize_smol_pos;
  initialization[KEY_COMMAND_OPTIONS] = "";
  initialization[KEY_EXPORT_ALL_ASCII] = true; // for testing, cellblender generates false as default
  initialization[KEY_MICROSCOPIC_REVERSIBILITY] = VALUE_OFF;
  initialization[KEY_TIME_STEP_MAX] = "";

  // reversed computation from mcell3's init_reactions
  pos_t vsd = sqrt_p(config.vacancy_search_dist2) * config.length_unit;
  initialization[KEY_VACANCY_SEARCH_DISTANCE] = f_to_str(vsd);

  initialization[KEY_SPACE_STEP] = "";
  initialization[KEY_SURFACE_GRID_DENSITY] = f_to_str(config.grid_density);

  // --- warnings ---
  Json::Value& warnings = initialization[KEY_WARNINGS];
  warnings[KEY_MISSED_REACTION_THRESHOLD] = "0.001";
  warnings[KEY_LIFETIME_TOO_SHORT] = VALUE_WARNING;
  warnings[KEY_LIFETIME_THRESHOLD] = "50";
  warnings[KEY_ALL_WARNINGS] = VALUE_INDIVIDUAL;
  warnings[KEY_MISSED_REACTIONS] = VALUE_WARNING;
  warnings[KEY_NEGATIVE_DIFFUSION_CONSTANT] = VALUE_WARNING;
  warnings[KEY_NEGATIVE_REACTION_RATE] = VALUE_WARNING;
  warnings[KEY_HIGH_PROBABILITY_THRESHOLD] = "1";
  warnings[KEY_DEGENERATE_POLYGONS] = VALUE_WARNING;
  warnings[KEY_USELESS_VOLUME_ORIENTATION] = VALUE_WARNING;

  warnings[KEY_HIGH_REACTION_PROBABILITY] =
      DMUtils::bool_to_warning_level(config.warnings.warn_on_bimol_rxn_probability_over_05_less_1);

  warnings[KEY_LARGE_MOLECULAR_DISPLACEMENT] = VALUE_WARNING;
  warnings[KEY_MISSING_SURFACE_ORIENTATION] = VALUE_ERROR;

  // --- notifications ---
  Json::Value& notifications = initialization[KEY_NOTIFICATIONS];

  notifications[KEY_SPECIES_REACTIONS_REPORT] = config.rxn_and_species_report;
  notifications[KEY_FILE_OUTPUT_REPORT] = false;
  notifications[KEY_ALL_NOTIFICATIONS] = VALUE_INDIVIDUAL;
  notifications[KEY_PROBABILITY_REPORT_THRESHOLD] = "0";
  notifications[KEY_BOX_TRIANGULATION_REPORT] = false;
  notifications[KEY_RELEASE_EVENT_REPORT] = true;
  notifications[KEY_PROGRESS_REPORT] = true;
  notifications[KEY_MOLECULE_COLLISION_REPORT] = false;
  notifications[KEY_ITERATION_REPORT] = true;
  notifications[KEY_FINAL_SUMMARY] = true;
  notifications[KEY_VARYING_PROBABILITY_REPORT] = config.notifications.rxn_probability_changed;
  notifications[KEY_PROBABILITY_REPORT] = VALUE_ON;
  notifications[KEY_PARTITION_LOCATION_REPORT] = false;
  notifications[KEY_DIFFUSION_CONSTANT_REPORT] = VALUE_BRIEF;

  // --- partitions ---
  Json::Value& partitions = initialization[KEY_PARTITIONS];
  // NOTE: mcell3_world_converter extends the partition info, so the result will be bigger than input
  // probably ok because we don't plan long-term MDL support without previous conversion
  partitions[KEY_INCLUDE] = true;
  partitions[KEY_RECURSION_FLAG] = false;

  const Vec3& origin = (config.partition0_llf * config.length_unit);
  pos_t length = config.partition_edge_length * config.length_unit;

  partitions[KEY_X_START] = f_to_str(origin.x);
  partitions[KEY_X_END] = f_to_str(origin.x + length);
  partitions[KEY_Y_START] = f_to_str(origin.y);
  partitions[KEY_Y_END] = f_to_str(origin.y + length);
  partitions[KEY_Z_START] = f_to_str(origin.z);
  partitions[KEY_Z_END] = f_to_str(origin.z + length);

  pos_t step = config.subpart_edge_length * config.length_unit;
  string step_str = f_to_str(step);
  partitions[KEY_X_STEP] = step_str;
  partitions[KEY_Y_STEP] = step_str;
  partitions[KEY_Z_STEP] = step_str;
}


void World::export_data_layout() const {

  Json::Value root;

  root[KEY_VERSION] = 2;
  root[KEY_MCELL4_MODE] = true;
  Json::Value& data_layout = root[KEY_DATA_LAYOUT];

  Json::Value dir;
  dir.append(VALUE_DIR);
  Json::Value dir_value;
  dir_value.append(".");
  dir.append(dir_value);
  data_layout.append(dir);

  // TODO: to go through viz and rxn outputs to
  // figure out what the directories should be, for now using the defaults
  Json::Value file_type;
  file_type.append(VALUE_FILE_TYPE);

  Json::Value file_type_contents;
  file_type_contents.append(VALUE_REACT_DATA);
  file_type_contents.append(VALUE_VIZ_DATA);

  file_type.append(file_type_contents);
  data_layout.append(file_type);

  Json::Value seed;
  seed.append(VALUE_SEED);
  Json::Value seed_value;
  seed_value.append(to_string(config.initial_seed));
  seed.append(seed_value);
  data_layout.append(seed);

  Json::StreamWriterBuilder wbuilder;
  wbuilder["indentation"] = " ";
  wbuilder.settings_["precision"] = 15; // this is the precision that is used by mdl_to_data_model.py script
  wbuilder.settings_["precisionType"] = "significant";
  std::string document = Json::writeString(wbuilder, root);

  // write result into a file
  ofstream res_file(DEFAULT_DATA_LAYOUT_FILENAME);
  if (res_file.is_open())
  {
    // maybe enable this message only in some verbose mode
    cout << "Generated file " << DEFAULT_DATA_LAYOUT_FILENAME << " for plotting in CellBlender.\n";
    res_file << document;
    res_file.close();
  }
  else {
    cout << "Unable to open file " << DEFAULT_DATA_LAYOUT_FILENAME << " for writing.\n";
  }
}

} // namespace mcell

