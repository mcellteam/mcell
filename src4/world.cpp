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

#include <fstream>

#include "rng.h" // MCell 3
#include "logging.h"

#include "world.h"
#include "viz_output_event.h"
#include "defragmentation_event.h"
#include "datamodel_defines.h"
#include "bng_data_to_datamodel_converter.h"


using namespace std;

const double USEC_IN_SEC = 1000000.0;

namespace MCell {

World::World()
  : bng_engine(config),
    total_iterations(0),
    seed_seq(0),
    next_wall_id(0),
    next_geometry_object_id(0),
    simulation_initialized(false),
    simulation_ended(false),
    previous_progress_report_time({0, 0}),
    previous_iteration(0),

    // temporary solution of callbacks
    wall_hit_callback(nullptr),
    wall_hit_callback_clientdata(nullptr),
    wall_hit_object_id(GEOMETRY_OBJECT_ID_INVALID)
{
  config.partition_edge_length = FLT_INVALID;
  config.subpartitions_per_partition_dimension = SUBPARTITIONS_PER_PARTITION_DIMENSION_DEFAULT;

  // although the same thing is called in init_simulation, not reseting it causes weird valdrind reports on
  // uninitialized variable
  reset_rusage(&sim_start_time);

#ifdef DEBUG_REACTIONS
  config.debug_reactions = true;
#endif
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



void World::init_counted_volumes() {
  assert(partitions.size() == 1);

  bool ok = CountedVolumesUtil::initialize_counted_volumes(this);
  if (!ok) {
    mcell_error("Processing of counted volumes failed, terminating.");
  }
}


void World::init_simulation() {
  assert(!simulation_initialized && "init_simulation must be called just once");

  if (get_all_species().get_count() == 0) {
    mcell_log("Error: there must be at lease one species!");
    exit(1);
  }

  // TODO: what do I need for initialization?
  config.init();
  stats.reset();

  init_fpu();

  init_counted_volumes();

  cout <<
      "Partition contains " <<  config.subpartitions_per_partition_dimension << "^3 subpartitions, " <<
      "subpartition size is " << config.subpartition_edge_length * config.length_unit << " microns.\n";
  assert(partitions.size() == 1 && "Initial partition must have been created, only 1 is allowed for now");

  // create defragmentation events
  DefragmentationEvent* defragmentation_event = new DefragmentationEvent(this);
  defragmentation_event->event_time = DEFRAGMENTATION_PERIODICITY;
  defragmentation_event->periodicity_interval = DEFRAGMENTATION_PERIODICITY;
  scheduler.schedule_event(defragmentation_event);

  // initialize timing
  previous_progress_report_time = {0, 0};

  rusage sim_start_time;
  reset_rusage(&sim_start_time);
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

  // flush and close count buffers
  for (CountBuffer& b: count_buffers) {
    b.flush_and_close();
  }

  // report final time
  rusage run_time;
  reset_rusage(&run_time);
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

  uint output_frequency = determine_output_frequency(total_iterations);

  // simulating 1000 iterations means to simulate iterations 0 .. 1000
  run_n_iterations(total_iterations + 1, output_frequency, true);

  end_simulation();
}


void World::dump() {
  stats.dump();

  get_all_species().dump(bng_engine.get_data());
  get_all_rxns().dump();

  // partitions
  for (Partition& p: partitions) {
    p.dump();
  }
}


void World::export_visualization_datamodel_to_dir(const char* prefix) const {
  // prefix should be the same directory that is used for viz_output,
  // e.g. ./viz_data/seed_0001/Scene

  stringstream path;
  path <<
      "4" << prefix << ".datamodel." <<
      VizOutputEvent::iterations_to_string(stats.get_current_iteration(), total_iterations) <<
      ".json";

  // create directories if needed
  ::make_parent_dir(path.str().c_str());

  export_visualization_datamodel(path.str().c_str());
}


void World::export_visualization_datamodel(const char* filename) const {

  Json::Value root;
  to_data_model(root);

  Json::StreamWriterBuilder wbuilder;
  wbuilder["indentation"] = " ";
  wbuilder.settings_["precision"] = 12;
  wbuilder.settings_["precisionType"] = "decimal";
  std::string document = Json::writeString(wbuilder, root);

  // write result into a file
  ofstream res_file(filename);
  if (res_file.is_open())
  {
    res_file << document;
    res_file.close();
  }
  else {
    cout << "Unable to open file " << filename << " for writing.\n";
  }
}


void World::to_data_model(Json::Value& root) const {
  Json::Value& mcell = root[KEY_MCELL];

  mcell[KEY_CELLBLENDER_VERSION] = VALUE_CELLBLENDER_VERSION;

  // generate geometry information
  bool first = true;
  for (const Partition& p: partitions) {
    p.to_data_model(mcell);
  }

  // generate species info
  BngDataToDatamodelConverter bng_converter;
  bng_converter.to_data_model(mcell, bng_engine);

  // add other default values, might need to generate this better
  Json::Value& materials = mcell[KEY_MATERIALS];
  Json::Value& material_dict = materials[KEY_MATERIAL_DICT];
  Json::Value& membrane = material_dict[KEY_VALUE_MEMBRANE];
  Json::Value& diffuse_color = membrane[KEY_DIFFUSE_COLOR];
  diffuse_color[KEY_R] = DEFAULT_OBJECT_COLOR_COMPONENT;
  diffuse_color[KEY_G] = DEFAULT_OBJECT_COLOR_COMPONENT;
  diffuse_color[KEY_B] = DEFAULT_OBJECT_COLOR_COMPONENT;
  diffuse_color[KEY_A] = DEFAULT_OBJECT_ALPHA;

  Json::Value& model_objects = mcell[KEY_MODEL_OBJECTS];
  DMUtil::json_add_version(model_objects, JSON_DM_VERSION_1330);
  Json::Value& model_object_list = model_objects[KEY_MODEL_OBJECT_LIST];

  // names of geometry objects need to be listed for the second time
  for (const Partition& p: partitions) {
    for (const GeometryObject& obj: p.get_geometry_objects()) {
      Json::Value model_object;
      model_object[KEY_PARENT_OBJECT] = "";
      model_object[KEY_DESCRIPTION] = "";
      model_object[KEY_OBJECT_SOURCE] = VALUE_BLENDER;
      model_object[KEY_DYNAMIC_DISPLAY_SOURCE] = "script";
      model_object[KEY_SCRIPT_NAME] = "";
      model_object[KEY_MEMBRANE_NAME] = "";
      model_object[KEY_DYNAMIC] = false;
      model_object[KEY_NAME] = DMUtil::remove_obj_name_prefix(obj.parent_name, obj.name);
      model_object_list.append(model_object);
    }
  }

  mcell[KEY_DEFINE_SURFACE_CLASSES] = Json::Value(Json::ValueType::objectValue); // empty dict

  Json::Value& mol_viz = mcell[KEY_MOL_VIZ];
  DMUtil::json_add_version(mol_viz, JSON_DM_VERSION_1700);
  mol_viz[KEY_MANUAL_SELECT_VIZ_DIR] = false;
  mol_viz[KEY_FILE_START_INDEX] = 0;
  mol_viz[KEY_SEED_LIST] = Json::Value(Json::arrayValue); // empty array

  Json::Value& color_list = mol_viz[KEY_COLOR_LIST];
  DMUtil::json_append_triplet(color_list, 0.8, 0.0, 0.0);
  DMUtil::json_append_triplet(color_list, 0.0, 0.8, 0.0);
  DMUtil::json_append_triplet(color_list, 0.0, 0.0, 0.8);
  DMUtil::json_append_triplet(color_list, 0.0, 0.8, 0.8);
  DMUtil::json_append_triplet(color_list, 0.8, 0.0, 0.8);
  DMUtil::json_append_triplet(color_list, 0.8, 0.8, 0.0);
  DMUtil::json_append_triplet(color_list, 1.0, 1.0, 1.0);
  DMUtil::json_append_triplet(color_list, 0.0, 0.0, 0.0);

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

  mcell[KEY_MODEL_LANGUAGE] = VALUE_MCELL3;

  mcell[KEY_PARAMETER_SYSTEM] = Json::Value(Json::ValueType::objectValue); // empty dict

  DMUtil::json_add_version(mcell, JSON_DM_VERSION_1300);

  mcell[KEY_INITIALIZATION] = Json::Value(Json::ValueType::objectValue); // empty dict

  Json::Value& blender_version = mcell[KEY_BLENDER_VERSION];
  blender_version.append(Json::Value(BLENDER_VERSION[0]));
  blender_version.append(Json::Value(BLENDER_VERSION[1]));
  blender_version.append(Json::Value(BLENDER_VERSION[2]));
}

} // namespace mcell

