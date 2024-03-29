/******************************************************************************
 *
 * Copyright (C) 2021 by
 * The Salk Institute for Biological Studies
 *
 * Use of this source code is governed by an MIT-style
 * license that can be found in the LICENSE file or at
 * https://opensource.org/licenses/MIT.
 *
******************************************************************************/

#include <sstream>
#include "api/pybind11_stl_include.h"
#include "api/python_export_utils.h"
#include "gen_config.h"
#include "api/api_config.h"
#include "api/rng_state.h"

namespace MCell {
namespace API {

void GenConfig::check_semantics() const {
}

void GenConfig::set_initialized() {
  if (is_set(initial_rng_state)) {
    initial_rng_state->set_initialized();
  }
  initialized = true;
}

void GenConfig::set_all_attributes_as_default_or_unset() {
  class_name = "Config";
  seed = 1;
  time_step = 1e-6;
  use_bng_units = false;
  surface_grid_density = 10000;
  interaction_radius = FLT_UNSET;
  intermembrane_interaction_radius = FLT_UNSET;
  vacancy_search_distance = 10;
  center_molecules_on_grid = false;
  partition_dimension = 10;
  initial_partition_origin = std::vector<double>();
  subpartition_dimension = 0.5;
  total_iterations = 1000000;
  check_overlapped_walls = true;
  reaction_class_cleanup_periodicity = 500;
  species_cleanup_periodicity = 10000;
  molecules_order_random_shuffle_periodicity = 10000;
  sort_molecules = false;
  memory_limit_gb = -1;
  initial_iteration = 0;
  initial_time = 0;
  initial_rng_state = nullptr;
  append_to_count_output_data = false;
  continue_after_sigalrm = false;
}

std::shared_ptr<Config> GenConfig::copy_config() const {
  std::shared_ptr<Config> res = std::make_shared<Config>(DefaultCtorArgType());
  res->class_name = class_name;
  res->seed = seed;
  res->time_step = time_step;
  res->use_bng_units = use_bng_units;
  res->surface_grid_density = surface_grid_density;
  res->interaction_radius = interaction_radius;
  res->intermembrane_interaction_radius = intermembrane_interaction_radius;
  res->vacancy_search_distance = vacancy_search_distance;
  res->center_molecules_on_grid = center_molecules_on_grid;
  res->partition_dimension = partition_dimension;
  res->initial_partition_origin = initial_partition_origin;
  res->subpartition_dimension = subpartition_dimension;
  res->total_iterations = total_iterations;
  res->check_overlapped_walls = check_overlapped_walls;
  res->reaction_class_cleanup_periodicity = reaction_class_cleanup_periodicity;
  res->species_cleanup_periodicity = species_cleanup_periodicity;
  res->molecules_order_random_shuffle_periodicity = molecules_order_random_shuffle_periodicity;
  res->sort_molecules = sort_molecules;
  res->memory_limit_gb = memory_limit_gb;
  res->initial_iteration = initial_iteration;
  res->initial_time = initial_time;
  res->initial_rng_state = initial_rng_state;
  res->append_to_count_output_data = append_to_count_output_data;
  res->continue_after_sigalrm = continue_after_sigalrm;

  return res;
}

std::shared_ptr<Config> GenConfig::deepcopy_config(py::dict) const {
  std::shared_ptr<Config> res = std::make_shared<Config>(DefaultCtorArgType());
  res->class_name = class_name;
  res->seed = seed;
  res->time_step = time_step;
  res->use_bng_units = use_bng_units;
  res->surface_grid_density = surface_grid_density;
  res->interaction_radius = interaction_radius;
  res->intermembrane_interaction_radius = intermembrane_interaction_radius;
  res->vacancy_search_distance = vacancy_search_distance;
  res->center_molecules_on_grid = center_molecules_on_grid;
  res->partition_dimension = partition_dimension;
  res->initial_partition_origin = initial_partition_origin;
  res->subpartition_dimension = subpartition_dimension;
  res->total_iterations = total_iterations;
  res->check_overlapped_walls = check_overlapped_walls;
  res->reaction_class_cleanup_periodicity = reaction_class_cleanup_periodicity;
  res->species_cleanup_periodicity = species_cleanup_periodicity;
  res->molecules_order_random_shuffle_periodicity = molecules_order_random_shuffle_periodicity;
  res->sort_molecules = sort_molecules;
  res->memory_limit_gb = memory_limit_gb;
  res->initial_iteration = initial_iteration;
  res->initial_time = initial_time;
  res->initial_rng_state = is_set(initial_rng_state) ? initial_rng_state->deepcopy_rng_state() : nullptr;
  res->append_to_count_output_data = append_to_count_output_data;
  res->continue_after_sigalrm = continue_after_sigalrm;

  return res;
}

bool GenConfig::__eq__(const Config& other) const {
  return
    seed == other.seed &&
    time_step == other.time_step &&
    use_bng_units == other.use_bng_units &&
    surface_grid_density == other.surface_grid_density &&
    interaction_radius == other.interaction_radius &&
    intermembrane_interaction_radius == other.intermembrane_interaction_radius &&
    vacancy_search_distance == other.vacancy_search_distance &&
    center_molecules_on_grid == other.center_molecules_on_grid &&
    partition_dimension == other.partition_dimension &&
    initial_partition_origin == other.initial_partition_origin &&
    subpartition_dimension == other.subpartition_dimension &&
    total_iterations == other.total_iterations &&
    check_overlapped_walls == other.check_overlapped_walls &&
    reaction_class_cleanup_periodicity == other.reaction_class_cleanup_periodicity &&
    species_cleanup_periodicity == other.species_cleanup_periodicity &&
    molecules_order_random_shuffle_periodicity == other.molecules_order_random_shuffle_periodicity &&
    sort_molecules == other.sort_molecules &&
    memory_limit_gb == other.memory_limit_gb &&
    initial_iteration == other.initial_iteration &&
    initial_time == other.initial_time &&
    (
      (is_set(initial_rng_state)) ?
        (is_set(other.initial_rng_state) ?
          (initial_rng_state->__eq__(*other.initial_rng_state)) : 
          false
        ) :
        (is_set(other.initial_rng_state) ?
          false :
          true
        )
     )  &&
    append_to_count_output_data == other.append_to_count_output_data &&
    continue_after_sigalrm == other.continue_after_sigalrm;
}

bool GenConfig::eq_nonarray_attributes(const Config& other, const bool ignore_name) const {
  return
    seed == other.seed &&
    time_step == other.time_step &&
    use_bng_units == other.use_bng_units &&
    surface_grid_density == other.surface_grid_density &&
    interaction_radius == other.interaction_radius &&
    intermembrane_interaction_radius == other.intermembrane_interaction_radius &&
    vacancy_search_distance == other.vacancy_search_distance &&
    center_molecules_on_grid == other.center_molecules_on_grid &&
    partition_dimension == other.partition_dimension &&
    true /*initial_partition_origin*/ &&
    subpartition_dimension == other.subpartition_dimension &&
    total_iterations == other.total_iterations &&
    check_overlapped_walls == other.check_overlapped_walls &&
    reaction_class_cleanup_periodicity == other.reaction_class_cleanup_periodicity &&
    species_cleanup_periodicity == other.species_cleanup_periodicity &&
    molecules_order_random_shuffle_periodicity == other.molecules_order_random_shuffle_periodicity &&
    sort_molecules == other.sort_molecules &&
    memory_limit_gb == other.memory_limit_gb &&
    initial_iteration == other.initial_iteration &&
    initial_time == other.initial_time &&
    (
      (is_set(initial_rng_state)) ?
        (is_set(other.initial_rng_state) ?
          (initial_rng_state->__eq__(*other.initial_rng_state)) : 
          false
        ) :
        (is_set(other.initial_rng_state) ?
          false :
          true
        )
     )  &&
    append_to_count_output_data == other.append_to_count_output_data &&
    continue_after_sigalrm == other.continue_after_sigalrm;
}

std::string GenConfig::to_str(const bool all_details, const std::string ind) const {
  std::stringstream ss;
  ss << get_object_name() << ": " <<
      "seed=" << seed << ", " <<
      "time_step=" << time_step << ", " <<
      "use_bng_units=" << use_bng_units << ", " <<
      "surface_grid_density=" << surface_grid_density << ", " <<
      "interaction_radius=" << interaction_radius << ", " <<
      "intermembrane_interaction_radius=" << intermembrane_interaction_radius << ", " <<
      "vacancy_search_distance=" << vacancy_search_distance << ", " <<
      "center_molecules_on_grid=" << center_molecules_on_grid << ", " <<
      "partition_dimension=" << partition_dimension << ", " <<
      "initial_partition_origin=" << vec_nonptr_to_str(initial_partition_origin, all_details, ind + "  ") << ", " <<
      "subpartition_dimension=" << subpartition_dimension << ", " <<
      "total_iterations=" << total_iterations << ", " <<
      "check_overlapped_walls=" << check_overlapped_walls << ", " <<
      "reaction_class_cleanup_periodicity=" << reaction_class_cleanup_periodicity << ", " <<
      "species_cleanup_periodicity=" << species_cleanup_periodicity << ", " <<
      "molecules_order_random_shuffle_periodicity=" << molecules_order_random_shuffle_periodicity << ", " <<
      "sort_molecules=" << sort_molecules << ", " <<
      "memory_limit_gb=" << memory_limit_gb << ", " <<
      "initial_iteration=" << initial_iteration << ", " <<
      "initial_time=" << initial_time << ", " <<
      "\n" << ind + "  " << "initial_rng_state=" << "(" << ((initial_rng_state != nullptr) ? initial_rng_state->to_str(all_details, ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "append_to_count_output_data=" << append_to_count_output_data << ", " <<
      "continue_after_sigalrm=" << continue_after_sigalrm;
  return ss.str();
}

py::class_<Config> define_pybinding_Config(py::module& m) {
  return py::class_<Config, std::shared_ptr<Config>>(m, "Config", "Class holds simulation configuration.")
      .def(
          py::init<
            const int,
            const double,
            const bool,
            const double,
            const double,
            const double,
            const double,
            const bool,
            const double,
            const std::vector<double>,
            const double,
            const double,
            const bool,
            const int,
            const int,
            const int,
            const bool,
            const int,
            const uint64_t,
            const double,
            std::shared_ptr<RngState>,
            const bool,
            const bool
          >(),
          py::arg("seed") = 1,
          py::arg("time_step") = 1e-6,
          py::arg("use_bng_units") = false,
          py::arg("surface_grid_density") = 10000,
          py::arg("interaction_radius") = FLT_UNSET,
          py::arg("intermembrane_interaction_radius") = FLT_UNSET,
          py::arg("vacancy_search_distance") = 10,
          py::arg("center_molecules_on_grid") = false,
          py::arg("partition_dimension") = 10,
          py::arg("initial_partition_origin") = std::vector<double>(),
          py::arg("subpartition_dimension") = 0.5,
          py::arg("total_iterations") = 1000000,
          py::arg("check_overlapped_walls") = true,
          py::arg("reaction_class_cleanup_periodicity") = 500,
          py::arg("species_cleanup_periodicity") = 10000,
          py::arg("molecules_order_random_shuffle_periodicity") = 10000,
          py::arg("sort_molecules") = false,
          py::arg("memory_limit_gb") = -1,
          py::arg("initial_iteration") = 0,
          py::arg("initial_time") = 0,
          py::arg("initial_rng_state") = nullptr,
          py::arg("append_to_count_output_data") = false,
          py::arg("continue_after_sigalrm") = false
      )
      .def("check_semantics", &Config::check_semantics)
      .def("__copy__", &Config::copy_config)
      .def("__deepcopy__", &Config::deepcopy_config, py::arg("memo"))
      .def("__str__", &Config::to_str, py::arg("all_details") = false, py::arg("ind") = std::string(""))
      .def("__eq__", &Config::__eq__, py::arg("other"))
      .def("dump", &Config::dump)
      .def_property("seed", &Config::get_seed, &Config::set_seed, "Random generator seed value.")
      .def_property("time_step", &Config::get_time_step, &Config::set_time_step, "Set the simulation time step to time_step seconds. 1e-6 (1us) is a common value. \nOne can set the time steps taken by individual molecules, but this \ntime step is still used as a default.\n")
      .def_property("use_bng_units", &Config::get_use_bng_units, &Config::set_use_bng_units, "When False (default), MCell uses traditional MCell units for bimolecular reaction rates are:\n * [M^-1*s^-1] for bimolecular reactions between either two volume molecules, a volume molecule \n               and a surface (molecule), \n * [um^2*N^-1*s^-1] bimolecular reactions between two surface molecules on the same surface.\nWhen True, BioNetGen units for bimolecular reaction rates are:\n * [um^3*N^-1*s^-1] for any bimolecular reactions. Surface-surface reaction rate conversion assumes 10nm membrane thickness\nBioNetGen units are compatible with BioNetGen's ODE, SSA, and PLA solvers given that seed species \nis copy number (N), these units are not compatible with NFSim. \nNo other units are affected by this setting.\n")
      .def_property("surface_grid_density", &Config::get_surface_grid_density, &Config::set_surface_grid_density, "Tile all surfaces so that they can hold molecules at N different positions per square micron.")
      .def_property("interaction_radius", &Config::get_interaction_radius, &Config::set_interaction_radius, "Diffusing volume molecules will interact with each other when\nthey get within N microns of each other. The default is\n1/sqrt(PI * Sigma_s) where Sigma_s is the surface grid density \n(default or user-specified).\n")
      .def_property("intermembrane_interaction_radius", &Config::get_intermembrane_interaction_radius, &Config::set_intermembrane_interaction_radius, "Diffusing surface molecules will interact with surface molecules on other\nwalls when they get within N microns of each other. The default is\n1/sqrt(PI * Sigma_s) where Sigma_s is the surface grid density \n(default or user-specified). \nWhen unset, the default value is computed as: \n1.0 / sqrt_f(MY_PI * surface_grid_density).\n")
      .def_property("vacancy_search_distance", &Config::get_vacancy_search_distance, &Config::set_vacancy_search_distance, "Rather internal, there is usually no need to change this value.\nUsed in dynamic geometry (see Model.apply_vertex_moves()). \nWhen a wall moves or its dimensions change, this is the maximum search distance \nuse when looking onto which tiles place the molecules on this wall. \nIf no empty tile is found within this distance, simulation fails.  \n    \n")
      .def_property("center_molecules_on_grid", &Config::get_center_molecules_on_grid, &Config::set_center_molecules_on_grid, "If set to True, then all molecules on a surface will be\nlocated exactly at the center of their grid element. If False, the\nmolecules will be randomly located when placed, and reactions\nwill take place at the location of the target (or the site of impact\nin the case of 3D molecule/surface reactions). \n")
      .def_property("partition_dimension", &Config::get_partition_dimension, &Config::set_partition_dimension, "All the simulated 3d space is placed in a partition. The partition is a cube and \nthis partition_dimension specifies the length of its edge in um.\n")
      .def_property("initial_partition_origin", &Config::get_initial_partition_origin, &Config::set_initial_partition_origin, py::return_value_policy::reference, "Optional placement of the initial partition in um, specifies the left, lower front \npoint. If not set, value -partition_dimension/2 is used for each of the dimensions \nplacing the center of the partition to (0, 0, 0).   \n")
      .def_property("subpartition_dimension", &Config::get_subpartition_dimension, &Config::set_subpartition_dimension, "Subpartition are spatial division of 3D space used to accelerate collision checking.\nIn general, partitions should be chosen to avoid having too many surfaces and molecules\nin one subpartition. \nIf there are few surfaces and/or molecules in a subvolume, it is advantageous to have the \nsubvolume as large as possible. Crossing partition boundaries takes a small amount of time, \nso it is rarely useful to have partitions more finely spaced than the average diffusion distance \nof the faster-moving molecules in the simulation.\n")
      .def_property("total_iterations", &Config::get_total_iterations, &Config::set_total_iterations, "Required for checkpointing so that the checkpointed model has information on\nthe intended total number of iterations. \nAlso used when generating visualization data files and also for other reporting uses. \nValue is truncated to an integer.\n")
      .def_property("check_overlapped_walls", &Config::get_check_overlapped_walls, &Config::set_check_overlapped_walls, "Enables check for overlapped walls. Overlapping walls can cause issues during \nsimulation such as a molecule escaping closed geometry when it hits two walls \nthat overlap. \n")
      .def_property("reaction_class_cleanup_periodicity", &Config::get_reaction_class_cleanup_periodicity, &Config::set_reaction_class_cleanup_periodicity, "Reaction class cleanup removes computed reaction classes for inactive species from memory.\nThis provides faster reaction lookup faster but when the same reaction class is \nneeded again, it must be recomputed.\n")
      .def_property("species_cleanup_periodicity", &Config::get_species_cleanup_periodicity, &Config::set_species_cleanup_periodicity, "Species cleanup removes inactive species from memory. It removes also all reaction classes \nthat reference it.\nThis provides faster addition of new species lookup faster but when the species is \nneeded again, it must be recomputed.\n")
      .def_property("molecules_order_random_shuffle_periodicity", &Config::get_molecules_order_random_shuffle_periodicity, &Config::set_molecules_order_random_shuffle_periodicity, "Randomly shuffle the order in which molecules are simulated.\nThis helps to overcome potential biases that may occur when \nmolecules are ordered e.g. by their species when simulation starts. \nThe first shuffling occurs at this iteration, i.e. no shuffle is done at iteration 0.\nSetting this parameter to 0 disables the shuffling.  \n")
      .def_property("sort_molecules", &Config::get_sort_molecules, &Config::set_sort_molecules, "Enables sorting of molecules for diffusion, this may improve cache locality and provide \nslightly better performance. \nProduces different results for the same seed when enabled because molecules are simulated \nin a different order. \n")
      .def_property("memory_limit_gb", &Config::get_memory_limit_gb, &Config::set_memory_limit_gb, "Sets memory limit in GB for simulation run. \nWhen this limit is hit, all buffers are flushed and simulation is terminated with an error.\n")
      .def_property("initial_iteration", &Config::get_initial_iteration, &Config::set_initial_iteration, "Initial iteration, used when resuming a checkpoint.")
      .def_property("initial_time", &Config::get_initial_time, &Config::set_initial_time, "Initial time in us, used when resuming a checkpoint.\nWill be truncated to be a multiple of time step.\n")
      .def_property("initial_rng_state", &Config::get_initial_rng_state, &Config::set_initial_rng_state, "Used for checkpointing, may contain state of the random number generator to be set \nafter initialization right before the first event is started. \nWhen not set, the set 'seed' value is used to initialize the random number generator.  \n")
      .def_property("append_to_count_output_data", &Config::get_append_to_count_output_data, &Config::set_append_to_count_output_data, "Used for checkpointing, instead of creating new files for Count observables data, \nnew values are appended to the existing files. If such files do not exist, new files are\ncreated.\n")
      .def_property("continue_after_sigalrm", &Config::get_continue_after_sigalrm, &Config::set_continue_after_sigalrm, "MCell registers a SIGALRM signal handler. When SIGALRM signal is received and \ncontinue_after_sigalrm is False, checkpoint is stored and simulation is terminated. \nWhen continue_after_sigalrm is True, checkpoint is stored and simulation continues.\nSIGALRM is not supported on Windows.\n")
    ;
}

std::string GenConfig::export_to_python(std::ostream& out, PythonExportContext& ctx) {
  if (!export_even_if_already_exported() && ctx.already_exported(this)) {
    return ctx.get_exported_name(this);
  }
  std::string exported_name = "config_" + std::to_string(ctx.postinc_counter("config"));
  if (!export_even_if_already_exported()) {
    ctx.add_exported(this, exported_name);
  }

  bool str_export = export_as_string_without_newlines();
  std::string nl = "";
  std::string ind = " ";
  std::stringstream ss;
  if (!str_export) {
    nl = "\n";
    ind = "    ";
    ss << exported_name << " = ";
  }
  ss << "m.Config(" << nl;
  if (seed != 1) {
    ss << ind << "seed = " << seed << "," << nl;
  }
  if (time_step != 1e-6) {
    ss << ind << "time_step = " << f_to_str(time_step) << "," << nl;
  }
  if (use_bng_units != false) {
    ss << ind << "use_bng_units = " << use_bng_units << "," << nl;
  }
  if (surface_grid_density != 10000) {
    ss << ind << "surface_grid_density = " << f_to_str(surface_grid_density) << "," << nl;
  }
  if (interaction_radius != FLT_UNSET) {
    ss << ind << "interaction_radius = " << f_to_str(interaction_radius) << "," << nl;
  }
  if (intermembrane_interaction_radius != FLT_UNSET) {
    ss << ind << "intermembrane_interaction_radius = " << f_to_str(intermembrane_interaction_radius) << "," << nl;
  }
  if (vacancy_search_distance != 10) {
    ss << ind << "vacancy_search_distance = " << f_to_str(vacancy_search_distance) << "," << nl;
  }
  if (center_molecules_on_grid != false) {
    ss << ind << "center_molecules_on_grid = " << center_molecules_on_grid << "," << nl;
  }
  if (partition_dimension != 10) {
    ss << ind << "partition_dimension = " << f_to_str(partition_dimension) << "," << nl;
  }
  if (initial_partition_origin != std::vector<double>() && !skip_vectors_export()) {
    ss << ind << "initial_partition_origin = " << export_vec_initial_partition_origin(out, ctx, exported_name) << "," << nl;
  }
  if (subpartition_dimension != 0.5) {
    ss << ind << "subpartition_dimension = " << f_to_str(subpartition_dimension) << "," << nl;
  }
  if (total_iterations != 1000000) {
    ss << ind << "total_iterations = " << f_to_str(total_iterations) << "," << nl;
  }
  if (check_overlapped_walls != true) {
    ss << ind << "check_overlapped_walls = " << check_overlapped_walls << "," << nl;
  }
  if (reaction_class_cleanup_periodicity != 500) {
    ss << ind << "reaction_class_cleanup_periodicity = " << reaction_class_cleanup_periodicity << "," << nl;
  }
  if (species_cleanup_periodicity != 10000) {
    ss << ind << "species_cleanup_periodicity = " << species_cleanup_periodicity << "," << nl;
  }
  if (molecules_order_random_shuffle_periodicity != 10000) {
    ss << ind << "molecules_order_random_shuffle_periodicity = " << molecules_order_random_shuffle_periodicity << "," << nl;
  }
  if (sort_molecules != false) {
    ss << ind << "sort_molecules = " << sort_molecules << "," << nl;
  }
  if (memory_limit_gb != -1) {
    ss << ind << "memory_limit_gb = " << memory_limit_gb << "," << nl;
  }
  if (initial_iteration != 0) {
    ss << ind << "initial_iteration = " << initial_iteration << "," << nl;
  }
  if (initial_time != 0) {
    ss << ind << "initial_time = " << f_to_str(initial_time) << "," << nl;
  }
  if (is_set(initial_rng_state)) {
    ss << ind << "initial_rng_state = " << initial_rng_state->export_to_python(out, ctx) << "," << nl;
  }
  if (append_to_count_output_data != false) {
    ss << ind << "append_to_count_output_data = " << append_to_count_output_data << "," << nl;
  }
  if (continue_after_sigalrm != false) {
    ss << ind << "continue_after_sigalrm = " << continue_after_sigalrm << "," << nl;
  }
  ss << ")" << nl << nl;
  if (!str_export) {
    out << ss.str();
    return exported_name;
  }
  else {
    return ss.str();
  }
}

std::string GenConfig::export_vec_initial_partition_origin(std::ostream& out, PythonExportContext& ctx, const std::string& parent_name) {
  // does not print the array itself to 'out' and returns the whole list
  std::stringstream ss;
  ss << "[";
  for (size_t i = 0; i < initial_partition_origin.size(); i++) {
    const auto& item = initial_partition_origin[i];
    if (i == 0) {
      ss << " ";
    }
    else if (i % 16 == 0) {
      ss << "\n  ";
    }
    ss << f_to_str(item) << ", ";
  }
  ss << "]";
  return ss.str();
}

} // namespace API
} // namespace MCell

