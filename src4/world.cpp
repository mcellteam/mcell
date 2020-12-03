/******************************************************************************
 *
 * Copyright (C) 2019, 2020 by
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
#include <run_n_iterations_end_event.h>
#include <sys/resource.h> // Linux include

#include <fstream>

#include "rng.h" // MCell 3
#include "logging.h"

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

#include "api/mol_wall_hit_info.h"
#include "api/geometry_object.h"
#include "api/model.h"

using namespace std;

const double USEC_IN_SEC = 1000000.0;

namespace MCell {

static double tousecs(timeval& t) {
  return (double)t.tv_sec * USEC_IN_SEC + (double)t.tv_usec;
}


static double tosecs(timeval& t) {
  return (double)t.tv_sec + (double)t.tv_usec/USEC_IN_SEC;
}


World::World(API::Callbacks& callbacks_)
  : bng_engine(config),
    callbacks(callbacks_),
    total_iterations(0),
    next_wall_id(0),
    next_region_id(0),
    next_geometry_object_id(0),
    simulation_initialized(false),
    simulation_ended(false),
    buffers_flushed(false),
    previous_progress_report_time({0, 0}),
    previous_iteration(0)
{
  config.partition_edge_length = FLT_INVALID;
  config.num_subpartitions_per_partition = SUBPARTITIONS_PER_PARTITION_DIMENSION_DEFAULT;

  // although the same thing is called in init_simulation, not reseting it causes weird valdrind reports on
  // uninitialized variable
  reset_rusage(&sim_start_time);
}


World::~World() {
  if (!buffers_flushed) {
    flush_buffers();
  }

  for (MolOrRxnCountEvent* e: unscheduled_count_events) {
    delete e;
  }
}

void World::init_fpu() {
  // empty
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


void World::recompute_species_flags() {
  // collect all MolOrRxnCountEvent to initialize species_flags_analyzer that is then used to set
  // certain species flags
  vector<BaseEvent*> count_events;
  scheduler.get_all_events_with_type_index(EVENT_TYPE_INDEX_MOL_OR_RXN_COUNT, count_events);
  species_flags_analyzer.initialize(count_events, unscheduled_count_events);
  get_all_species().recompute_species_flags(get_all_rxns(), &species_flags_analyzer);
}

void World::init_counted_volumes() {
  assert(partitions.size() == 1);

  bool ok = CountedVolumesUtil::initialize_counted_volumes(this, config.has_intersecting_counted_objects);
  if (!ok) {
    mcell_error("Processing of counted volumes failed, terminating.");
  }

  partitions[PARTITION_ID_INITIAL].initialize_all_waypoints();
}


void World::init_simulation() {

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

  recompute_species_flags();

  stats.reset();

  init_fpu();

  init_counted_volumes();

  cout <<
      "Partition contains " <<  config.num_subpartitions_per_partition << "^3 subpartitions, " <<
      "subpartition size is " << config.subpartition_edge_length * config.length_unit << " microns.\n";
  assert(partitions.size() == 1 && "Initial partition must have been created, only 1 is allowed for now");

  // create event that diffuses molecules
  DiffuseReactEvent* event = new DiffuseReactEvent(this);
  event->event_time = TIME_SIMULATION_START;
  scheduler.schedule_event(event);

  // create defragmentation events
  DefragmentationEvent* defragmentation_event = new DefragmentationEvent(this);
  defragmentation_event->event_time = DEFRAGMENTATION_PERIODICITY;
  defragmentation_event->periodicity_interval = DEFRAGMENTATION_PERIODICITY;
  scheduler.schedule_event(defragmentation_event);

  // create rxn class cleanup events
  RxnClassCleanupEvent* rxn_class_cleanup_event = new RxnClassCleanupEvent(this);
  rxn_class_cleanup_event->event_time = RXN_CLASS_CLEANUP_PERIODICITY;
  rxn_class_cleanup_event->periodicity_interval = RXN_CLASS_CLEANUP_PERIODICITY;
  scheduler.schedule_event(rxn_class_cleanup_event);

  SpeciesCleanupEvent* species_cleanup_event = new SpeciesCleanupEvent(this);
  species_cleanup_event->event_time = SPECIES_CLEANUP_PERIODICITY;
  species_cleanup_event->periodicity_interval = SPECIES_CLEANUP_PERIODICITY;
  scheduler.schedule_event(species_cleanup_event);

  // create subpart sorting events
  if (config.sort_mols_by_subpart) {
    SortMolsBySubpartEvent* sort_event = new SortMolsBySubpartEvent(this);
    sort_event->event_time = 0;
    sort_event->periodicity_interval = SORT_MOLS_BY_SUBPART_PERIODICITY;
    scheduler.schedule_event(sort_event);
  }

  // initialize timing
  previous_progress_report_time = {0, 0};

  rusage sim_start_time;
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

  simulation_initialized = true;
}


void World::run_n_iterations(const uint64_t num_iterations, const uint64_t output_frequency, const bool terminate_last_iteration_after_viz_output) {

  if (!simulation_initialized) {
    init_simulation();
  }

  uint64_t& current_iteration = stats.get_current_iteration();

  // create events that are used to check whether simulation should end, also serves as a barrier
  // also serves as a simulation barrier to not to do diffusion after this point in time
  RunNIterationsEndEvent* run_n_iterations_end_event = new RunNIterationsEndEvent();
  run_n_iterations_end_event->event_time = current_iteration + num_iterations;
  run_n_iterations_end_event->periodicity_interval = 1; // these markers are inserted into every time step
  scheduler.schedule_event(run_n_iterations_end_event);

  if (current_iteration == 0) {
    cout << "Iterations: " << current_iteration << " of " << total_iterations << "\n";
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
        cout << "Iterations: " << current_iteration << " of " << total_iterations;

        timeval current_progress_report_time;
        gettimeofday(&current_progress_report_time, NULL);
        if (previous_progress_report_time.tv_usec > 0) {
          double time_diff = tousecs(current_progress_report_time) - tousecs(previous_progress_report_time);
          time_diff /= (double)output_frequency;
          cout << " (" << 1000000.0/time_diff << " iter/sec)";
        }
        previous_progress_report_time = current_progress_report_time;

        cout << " " << bng_engine.get_stats_report();

        cout << "\n";
        cout.flush(); // flush is required so that CellBlender can display progress
      }

      previous_iteration = current_iteration;
    }

#ifdef DEBUG_SCHEDULER
    cout << "After it: " << current_iteration << ", time: " << time << "\n";
#endif

    // also terminate if this was the last iteration and we hit an event that represents a check for the
    // end of the simulation
    if (terminate_last_iteration_after_viz_output &&
       event_info.type_index == EVENT_TYPE_INDEX_SIMULATION_END_CHECK
    ) {
      assert(current_iteration == this_run_first_iteration + num_iterations - 1);
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


void World::flush_buffers() {
  assert(!buffers_flushed && "Buffers can be flushed only once");

  // flush and close count buffers
  for (CountBuffer& b: count_buffers) {
    b.flush_and_close();
  }
  buffers_flushed = true;
}


void World::end_simulation(const bool run_up_to_last_count_and_viz_count_events, const bool print_final_report) {
  if (simulation_ended) {
    // already called, do nothing
    return;
  }

  if (run_up_to_last_count_and_viz_count_events) {
    // executes all
    run_n_iterations(1, determine_output_frequency(total_iterations), true);
  }

  flush_buffers();

  if (print_final_report) {
    cout << "Iteration " << stats.get_current_iteration() << ", simulation finished successfully\n";

    stats.dump();

    // report final time
    rusage run_time;
    reset_rusage(&run_time);
    getrusage(RUSAGE_SELF, &run_time);
    cout << "Simulation CPU time = "
      << tosecs(run_time.ru_utime) - tosecs(sim_start_time.ru_utime) <<  "(user) and "
      << tosecs(run_time.ru_stime) - tosecs(sim_start_time.ru_stime) <<  "(system)\n";
  }

  simulation_ended = true;
}


void World::run_simulation(const bool dump_initial_state) {

  // do initialization, also insert
  // defragmentation and end simulation event
  init_simulation();

  if (dump_initial_state) {
    dump();
  }

  uint output_frequency = World::determine_output_frequency(total_iterations);

  // simulating 1000 iterations means to simulate iterations 0 .. 1000
  run_n_iterations(total_iterations + 1, output_frequency, true);

  end_simulation();
}


void World::dump() {
  config.dump();
  stats.dump();

  bng_engine.get_data().dump();
  get_all_species().dump();
  get_all_rxns().dump(true);


  // partitions
  for (Partition& p: partitions) {
    p.dump();
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

  set<species_id_t> affected_species;
  const auto& users = rxn->get_rxn_classed_where_used();
  for (const auto& u: users) {
    assert(u->specific_reactants.size() == 1);
    affected_species.insert(u->specific_reactants[0].species_id);
  }

  // and then reset unimol time for each molecule of that species
  for (Partition& p: partitions) {
    for (Molecule& m: p.get_molecules()) {
      if (affected_species.count(m.species_id) != 0) {
        // new unimol time will be computed when the molecule is diffused
        // the next time (we cannot change it right away for molecules that
        // have longer timestep anyway because they were already diffused to the future)
        m.unimol_rx_time = TIME_INVALID;
        m.clear_flag(MOLECULE_FLAG_RESCHEDULE_UNIMOL_RXN_ON_NEXT_RXN_RATE_UPDATE);
        m.set_flag(MOLECULE_FLAG_SCHEDULE_UNIMOL_RXN);
      }
    }
  }
}



std::string World::export_releases_to_bngl_seed_species(
    std::ostream& parameters, std::ostream& seed_species) const {
  seed_species << BNG::BEGIN_SEED_SPECIES << "\n";

  parameters << "\n" << BNG::IND << "# seed species counts\n";

  vector<const BaseEvent*> release_events;
  scheduler.get_all_events_with_type_index(EVENT_TYPE_INDEX_RELEASE, release_events);

  for (size_t i = 0; i < release_events.size(); i++) {
    const ReleaseEvent* re = dynamic_cast<const ReleaseEvent*>(release_events[i]);
    assert(re != nullptr);

    if (re->release_shape == ReleaseShape::INITIAL_SURF_REGION) {
      // TODO: check whether there are initial releases, ignoring it for now because it is not often used
      continue;
    }

    const string & err_suffix = ", error for " + re->release_site_name + ".";
    if (re->release_shape != ReleaseShape::REGION) {
      return "Only region release shapes are currently supported for BNGL export" + err_suffix;
    }
    if (re->release_number_method != ReleaseNumberMethod::ConstNum) {
      return "Only constant release number releases are currently supported for BNGL export" + err_suffix;
    }
    if (re->region_expr_root->op != RegionExprOperator::Leaf) {
      return "Only simple release regions are currently supported for BNGL export" + err_suffix;
    }
    if (re->event_time != 0) {
      return "Only releases for time 0 are currently supported for BNGL export" + err_suffix;
    }
    if (re->orientation != ORIENTATION_NONE) {
      return "Only releases for volume molecules are currently supported for BNGL export" + err_suffix;
    }
    if (re->needs_release_pattern()) {
      return "Releases with release patterns are not currently supported for BNGL export" + err_suffix;
    }

    // create parameter for the count
    string seed_count_name = "seed_count_" + to_string(i);
    parameters << BNG::IND << seed_count_name << " " << to_string(re->release_number) << "\n";

    // and line in seed species, for now whole objects are representing compartments
    const Region& region = get_region(re->region_expr_root->region_id);
    if (DMUtil::get_region_name(region.name) != "ALL") {
      return "Compartments that do not span the whole object are not supported yet" + err_suffix;
    }
    const GeometryObject& obj = get_geometry_object(region.geometry_object_id);
    const BNG::Species& species = bng_engine.get_all_species().get(re->species_id);
    seed_species << BNG::IND <<
        "@" <<  obj.name << ":" << species.name << " " << seed_count_name << "\n";
  }

  seed_species << BNG::END_SEED_SPECIES << "\n";
  return "";
}


std::string World::export_counts_to_bngl_observables(std::ostream& observables) const {
  observables << BNG::BEGIN_OBSERVABLES << "\n";

  vector<const BaseEvent*> count_events;
  scheduler.get_all_events_with_type_index(EVENT_TYPE_INDEX_MOL_OR_RXN_COUNT, count_events);

  for (size_t i = 0; i < count_events.size(); i++) {
    const MolOrRxnCountEvent* ce = dynamic_cast<const MolOrRxnCountEvent*>(count_events[i]);
    assert(ce != nullptr);

    for (const MolOrRxnCountItem& item: ce->mol_count_items) {
      // get observable name from filename
      const string& path = get_count_buffer(item.buffer_id).get_filename();
      size_t slash_pos = path.find_last_of("/\\");
      release_assert(slash_pos != string::npos);
      size_t dot_pos = path.rfind('.');
      release_assert(dot_pos != string::npos);
      string name = path.substr(slash_pos + 1, dot_pos - slash_pos - 1);

      const string& err_suffix = ", error for " + name + ".";

      if (item.terms.size() != 1 || item.multiplier != 1.0) {
        return "Observable expressions other than a single pattern are not yet supported by BNGL export" + err_suffix;
      }
      const MolOrRxnCountTerm& term = item.terms[0];
      if (term.is_rxn_count()) {
        return "Reaction counts are not supported by BNGL export" + err_suffix;
      }
      if (term.type == CountType::PresentOnSurfaceRegion) {
        return "Surface molecule counts are yet not supported by BNGL export" + err_suffix;
      }

      string compartment_prefix = "";
      if (term.type == CountType::EnclosedInVolumeRegion) {
        compartment_prefix = "@" + get_geometry_object(term.geometry_object_id).name + ":";
      }

      // Species or Molecules
      string type;
      if (term.species_pattern_type == SpeciesPatternType::SpeciesPattern) {
        type = BNG::OBSERVABLE_SPECIES;
      }
      else if (term.species_pattern_type == SpeciesPatternType::MoleculesPattern) {
        type = BNG::OBSERVABLE_MOLECULES;
      }
      else {
        release_assert("SpeciesId type should not be used here.");
      }

      string pattern = term.species_molecules_pattern.to_str(false);

      observables << BNG::IND <<
          type << " " << name << " " << compartment_prefix << pattern << " " << "\n";
    }
  }

  observables << BNG::END_OBSERVABLES << "\n";
  return "";
}


// returns empty string if everything went well, nonempty string with error message
std::string World::export_to_bngl(const std::string& file_name) const {

  ofstream out;
  out.open(file_name);
  if (!out.is_open()) {
    return "Could not open output file " + file_name + ".";
  }


  const Partition& p = get_partition(PARTITION_ID_INITIAL);
  if (p.get_geometry_objects().size() != 1) {
    return "The only supported models currently are those that have exactly 1 geometry object.";
  }
  const GeometryObject& obj = p.get_geometry_objects()[0];

  float_t volume_internal_units = CountedVolumesUtil::get_geometry_object_volume(this, obj);
  if (volume_internal_units == FLT_INVALID) {
    return "Compartment object " + obj.name + " is not watertight and its volume cannot be computed.";
  }

  float_t volume = volume_internal_units * pow(config.length_unit, 3);

  stringstream parameters;
  stringstream molecule_types;
  stringstream reaction_rules;
  parameters << BNG::BEGIN_PARAMETERS << "\n";

  parameters << BNG::IND << "# general parameters\n";
  parameters << BNG::IND << BNG::ITERATIONS << " " << total_iterations << "\n";
  parameters << BNG::IND << BNG::MCELL_TIME_STEP << " " << f_to_str(config.time_unit) << "\n";

  string err_msg = bng_engine.export_to_bngl(
      parameters, molecule_types, reaction_rules, volume);
  if (err_msg != "") {
    out.close();
    return err_msg;
  }

  // seed species
  stringstream seed_species;
  err_msg = export_releases_to_bngl_seed_species(parameters, seed_species);
  if (err_msg != "") {
    out.close();
    return err_msg;
  }

  stringstream observables;
  err_msg = export_counts_to_bngl_observables(observables);
  if (err_msg != "") {
    out.close();
    return err_msg;
  }

  parameters << BNG::END_PARAMETERS << "\n";

  out << parameters.str() << "\n";
  out << molecule_types.str() << "\n";

  // compartments
  out << BNG::BEGIN_COMPARTMENTS << "\n";
  out << BNG::IND << obj.name << " 3 " << BNG::PARAM_V << " * 1e15 # volume in fL (um^3)\n";
  out << BNG::END_COMPARTMENTS << "\n\n";

  out << seed_species.str() << "\n";
  out << observables.str() << "\n";

  out << reaction_rules.str() << "\n";

  // NFSim/ODE
  // simulate({method=>"nf",seed=>1,gml=>1000000,t_end=>1e-3,n_steps=>1000})

  out.close();

  return "";
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
  ::make_parent_dir(path.str().c_str());

  export_data_model(path.str().c_str(), only_for_viz);
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
  DMUtil::add_version(mcell, VER_DM_2017_06_23_1300);

  initialization_to_data_model(mcell);

  // generate geometry information

  // first create empty model_objects section (may be filled-in by
  // GeometryObject::to_data_model_as_model_object
  Json::Value& model_objects = mcell[KEY_MODEL_OBJECTS];
  DMUtil::add_version(model_objects, VER_DM_2018_01_11_1330);
  Json::Value& model_object_list = model_objects[KEY_MODEL_OBJECT_LIST];
  model_object_list = Json::Value(Json::arrayValue);

  // then dump all partition data
  bool first = true;
  for (const Partition& p: partitions) {
    p.to_data_model(mcell);
  }

  if (!only_for_viz) {
    // base information for reaction_data_output must be set even when there are no such events
    Json::Value& reaction_data_output = mcell[KEY_REACTION_DATA_OUTPUT];
    DMUtil::add_version(reaction_data_output, VER_DM_2016_03_15_1800);
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
  Json::Value& membrane = material_dict[KEY_VALUE_MEMBRANE];
  Json::Value& diffuse_color = membrane[KEY_DIFFUSE_COLOR];
  diffuse_color[KEY_R] = DEFAULT_OBJECT_COLOR_COMPONENT;
  diffuse_color[KEY_G] = DEFAULT_OBJECT_COLOR_COMPONENT;
  diffuse_color[KEY_B] = DEFAULT_OBJECT_COLOR_COMPONENT;
  diffuse_color[KEY_A] = DEFAULT_OBJECT_ALPHA;

  // diverse settings not read from the mcell4 state

  // --- mol_viz ---
  Json::Value& mol_viz = mcell[KEY_MOL_VIZ];
  DMUtil::add_version(mol_viz, VER_DM_2015_04_13_1700);
  mol_viz[KEY_MANUAL_SELECT_VIZ_DIR] = false;
  mol_viz[KEY_FILE_START_INDEX] = 0;
  mol_viz[KEY_SEED_LIST] = Json::Value(Json::arrayValue); // empty array

  Json::Value& color_list = mol_viz[KEY_COLOR_LIST];
  DMUtil::append_triplet(color_list, 0.8, 0.0, 0.0);
  DMUtil::append_triplet(color_list, 0.0, 0.8, 0.0);
  DMUtil::append_triplet(color_list, 0.0, 0.0, 0.8);
  DMUtil::append_triplet(color_list, 0.0, 0.8, 0.8);
  DMUtil::append_triplet(color_list, 0.8, 0.0, 0.8);
  DMUtil::append_triplet(color_list, 0.8, 0.8, 0.0);
  DMUtil::append_triplet(color_list, 1.0, 1.0, 1.0);
  DMUtil::append_triplet(color_list, 0.0, 0.0, 0.0);

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
  DMUtil::add_version(scripting, VER_DM_2017_11_30_1830);
  scripting[KEY_SCRIPTING_LIST] = Json::Value(Json::arrayValue);
  scripting[KEY_SCRIPT_TEXTS] = Json::Value(Json::objectValue);
  scripting[KEY_DM_INTERNAL_FILE_NAME] = "";
  scripting[KEY_FORCE_PROPERTY_UPDATE] = true;
  scripting[KEY_DM_EXTERNAL_FILE_NAME] = "";
  scripting[KEY_IGNORE_CELLBLENDER_DATA] = false;

  // --- simulation_control ---
  Json::Value& simulation_control = mcell[KEY_SIMULATION_CONTROL];
  simulation_control[KEY_EXPORT_FORMAT] = VALUE_MCELL_MDL_MODULAR;

  mcell[KEY_MODEL_LANGUAGE] = VALUE_MCELL3;

  Json::Value& blender_version = mcell[KEY_BLENDER_VERSION];
  blender_version.append(Json::Value(BLENDER_VERSION[0]));
  blender_version.append(Json::Value(BLENDER_VERSION[1]));
  blender_version.append(Json::Value(BLENDER_VERSION[2]));
}


void World::initialization_to_data_model(Json::Value& mcell_node) const {
  // only setting defaults for now, most of these values are not used in mcell4

  // --- initialization ---
  Json::Value& initialization = mcell_node[KEY_INITIALIZATION];
  DMUtil::add_version(initialization, VER_DM_2017_11_18_0130);

  // time step will most probably use rounded values, therefore we don't have to use full precision here
  initialization[KEY_TIME_STEP] = f_to_str(config.time_unit, 8);
  initialization[KEY_ITERATIONS] = to_string(total_iterations);

  initialization[KEY_INTERACTION_RADIUS] = "";
  initialization[KEY_ACCURATE_3D_REACTIONS] = true;
  initialization[KEY_RADIAL_SUBDIVISIONS] = "";
  initialization[KEY_RADIAL_DIRECTIONS] = "";
  initialization[KEY_CENTER_MOLECULES_ON_GRID] = !config.randomize_smol_pos;
  initialization[KEY_COMMAND_OPTIONS] = "";
  initialization[KEY_EXPORT_ALL_ASCII] = true; // for testing, cellblender generates false as default
  initialization[KEY_MICROSCOPIC_REVERSIBILITY] = VALUE_OFF;
  initialization[KEY_TIME_STEP_MAX] = "";

  // reversed computation from mcell3's init_reactions
  float_t vsd = sqrt(config.vacancy_search_dist2) * config.length_unit;
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
  warnings[KEY_HIGH_REACTION_PROBABILITY] = VALUE_IGNORED;
  warnings[KEY_LARGE_MOLECULAR_DISPLACEMENT] = VALUE_WARNING;
  warnings[KEY_MISSING_SURFACE_ORIENTATION] = VALUE_ERROR;

  // --- notifications ---
  Json::Value& notifications = initialization[KEY_NOTIFICATIONS];
  notifications[KEY_FILE_OUTPUT_REPORT] = false;
  notifications[KEY_ALL_NOTIFICATIONS] = VALUE_INDIVIDUAL;
  notifications[KEY_PROBABILITY_REPORT_THRESHOLD] = "0";
  notifications[KEY_BOX_TRIANGULATION_REPORT] = false;
  notifications[KEY_RELEASE_EVENT_REPORT] = true;
  notifications[KEY_PROGRESS_REPORT] = true;
  notifications[KEY_MOLECULE_COLLISION_REPORT] = false;
  notifications[KEY_ITERATION_REPORT] = true;
  notifications[KEY_FINAL_SUMMARY] = true;
  notifications[KEY_VARYING_PROBABILITY_REPORT] = true;
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
  float_t length = config.partition_edge_length * config.length_unit;

  partitions[KEY_X_START] = f_to_str(origin.x);
  partitions[KEY_X_END] = f_to_str(origin.x + length);
  partitions[KEY_Y_START] = f_to_str(origin.y);
  partitions[KEY_Y_END] = f_to_str(origin.y + length);
  partitions[KEY_Z_START] = f_to_str(origin.z);
  partitions[KEY_Z_END] = f_to_str(origin.z + length);

  float_t step = config.subpartition_edge_length * config.length_unit;
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

