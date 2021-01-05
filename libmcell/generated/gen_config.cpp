/******************************************************************************
 *
 * Copyright (C) 2020 by
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

#include <sstream>
#include "libs/pybind11/include/pybind11/stl.h"
#include "api/python_export_utils.h"
#include "gen_config.h"
#include "api/config.h"
#include "api/rng_state.h"

namespace MCell {
namespace API {

void GenConfig::check_semantics() const {
}

void GenConfig::set_initialized() {
  if (is_set(rng_state)) {
    rng_state->set_initialized();
  }
  initialized = true;
}

void GenConfig::set_all_attributes_as_default_or_unset() {
  class_name = "Config";
  seed = 1;
  time_step = 1e-6;
  surface_grid_density = 10000;
  interaction_radius = FLT_UNSET;
  intermembrane_interaction_radius = FLT_UNSET;
  vacancy_search_distance = 10;
  center_molecules_on_grid = false;
  initial_partition_origin = std::vector<float_t>();
  partition_dimension = 10;
  subpartition_dimension = 0.5;
  total_iterations = 1000000;
  check_overlapped_walls = true;
  reaction_class_cleanup_periodicity = 500;
  species_cleanup_periodicity = 10000;
  sort_molecules = false;
  memory_limit_gb = -1;
  initial_iteration = 0;
  initial_time = 0;
  rng_state = nullptr;
}

bool GenConfig::__eq__(const Config& other) const {
  return
    seed == other.seed &&
    time_step == other.time_step &&
    surface_grid_density == other.surface_grid_density &&
    interaction_radius == other.interaction_radius &&
    intermembrane_interaction_radius == other.intermembrane_interaction_radius &&
    vacancy_search_distance == other.vacancy_search_distance &&
    center_molecules_on_grid == other.center_molecules_on_grid &&
    initial_partition_origin == other.initial_partition_origin &&
    partition_dimension == other.partition_dimension &&
    subpartition_dimension == other.subpartition_dimension &&
    total_iterations == other.total_iterations &&
    check_overlapped_walls == other.check_overlapped_walls &&
    reaction_class_cleanup_periodicity == other.reaction_class_cleanup_periodicity &&
    species_cleanup_periodicity == other.species_cleanup_periodicity &&
    sort_molecules == other.sort_molecules &&
    memory_limit_gb == other.memory_limit_gb &&
    initial_iteration == other.initial_iteration &&
    initial_time == other.initial_time &&
    (
      (is_set(rng_state)) ?
        (is_set(other.rng_state) ?
          (rng_state->__eq__(*other.rng_state)) : 
          false
        ) :
        (is_set(other.rng_state) ?
          false :
          true
        )
     ) ;
}

bool GenConfig::eq_nonarray_attributes(const Config& other, const bool ignore_name) const {
  return
    seed == other.seed &&
    time_step == other.time_step &&
    surface_grid_density == other.surface_grid_density &&
    interaction_radius == other.interaction_radius &&
    intermembrane_interaction_radius == other.intermembrane_interaction_radius &&
    vacancy_search_distance == other.vacancy_search_distance &&
    center_molecules_on_grid == other.center_molecules_on_grid &&
    true /*initial_partition_origin*/ &&
    partition_dimension == other.partition_dimension &&
    subpartition_dimension == other.subpartition_dimension &&
    total_iterations == other.total_iterations &&
    check_overlapped_walls == other.check_overlapped_walls &&
    reaction_class_cleanup_periodicity == other.reaction_class_cleanup_periodicity &&
    species_cleanup_periodicity == other.species_cleanup_periodicity &&
    sort_molecules == other.sort_molecules &&
    memory_limit_gb == other.memory_limit_gb &&
    initial_iteration == other.initial_iteration &&
    initial_time == other.initial_time &&
    (
      (is_set(rng_state)) ?
        (is_set(other.rng_state) ?
          (rng_state->__eq__(*other.rng_state)) : 
          false
        ) :
        (is_set(other.rng_state) ?
          false :
          true
        )
     ) ;
}

std::string GenConfig::to_str(const std::string ind) const {
  std::stringstream ss;
  ss << get_object_name() << ": " <<
      "seed=" << seed << ", " <<
      "time_step=" << time_step << ", " <<
      "surface_grid_density=" << surface_grid_density << ", " <<
      "interaction_radius=" << interaction_radius << ", " <<
      "intermembrane_interaction_radius=" << intermembrane_interaction_radius << ", " <<
      "vacancy_search_distance=" << vacancy_search_distance << ", " <<
      "center_molecules_on_grid=" << center_molecules_on_grid << ", " <<
      "initial_partition_origin=" << vec_nonptr_to_str(initial_partition_origin, ind + "  ") << ", " <<
      "partition_dimension=" << partition_dimension << ", " <<
      "subpartition_dimension=" << subpartition_dimension << ", " <<
      "total_iterations=" << total_iterations << ", " <<
      "check_overlapped_walls=" << check_overlapped_walls << ", " <<
      "reaction_class_cleanup_periodicity=" << reaction_class_cleanup_periodicity << ", " <<
      "species_cleanup_periodicity=" << species_cleanup_periodicity << ", " <<
      "sort_molecules=" << sort_molecules << ", " <<
      "memory_limit_gb=" << memory_limit_gb << ", " <<
      "initial_iteration=" << initial_iteration << ", " <<
      "initial_time=" << initial_time << ", " <<
      "\n" << ind + "  " << "rng_state=" << "(" << ((rng_state != nullptr) ? rng_state->to_str(ind + "  ") : "null" ) << ")";
  return ss.str();
}

py::class_<Config> define_pybinding_Config(py::module& m) {
  return py::class_<Config, std::shared_ptr<Config>>(m, "Config")
      .def(
          py::init<
            const int,
            const float_t,
            const float_t,
            const float_t,
            const float_t,
            const float_t,
            const bool,
            const std::vector<float_t>,
            const float_t,
            const float_t,
            const float_t,
            const bool,
            const int,
            const int,
            const bool,
            const int,
            const uint64_t,
            const float_t,
            std::shared_ptr<RngState>
          >(),
          py::arg("seed") = 1,
          py::arg("time_step") = 1e-6,
          py::arg("surface_grid_density") = 10000,
          py::arg("interaction_radius") = FLT_UNSET,
          py::arg("intermembrane_interaction_radius") = FLT_UNSET,
          py::arg("vacancy_search_distance") = 10,
          py::arg("center_molecules_on_grid") = false,
          py::arg("initial_partition_origin") = std::vector<float_t>(),
          py::arg("partition_dimension") = 10,
          py::arg("subpartition_dimension") = 0.5,
          py::arg("total_iterations") = 1000000,
          py::arg("check_overlapped_walls") = true,
          py::arg("reaction_class_cleanup_periodicity") = 500,
          py::arg("species_cleanup_periodicity") = 10000,
          py::arg("sort_molecules") = false,
          py::arg("memory_limit_gb") = -1,
          py::arg("initial_iteration") = 0,
          py::arg("initial_time") = 0,
          py::arg("rng_state") = nullptr
      )
      .def("check_semantics", &Config::check_semantics)
      .def("__str__", &Config::to_str, py::arg("ind") = std::string(""))
      .def("__eq__", &Config::__eq__, py::arg("other"))
      .def("dump", &Config::dump)
      .def_property("seed", &Config::get_seed, &Config::set_seed)
      .def_property("time_step", &Config::get_time_step, &Config::set_time_step)
      .def_property("surface_grid_density", &Config::get_surface_grid_density, &Config::set_surface_grid_density)
      .def_property("interaction_radius", &Config::get_interaction_radius, &Config::set_interaction_radius)
      .def_property("intermembrane_interaction_radius", &Config::get_intermembrane_interaction_radius, &Config::set_intermembrane_interaction_radius)
      .def_property("vacancy_search_distance", &Config::get_vacancy_search_distance, &Config::set_vacancy_search_distance)
      .def_property("center_molecules_on_grid", &Config::get_center_molecules_on_grid, &Config::set_center_molecules_on_grid)
      .def_property("initial_partition_origin", &Config::get_initial_partition_origin, &Config::set_initial_partition_origin)
      .def_property("partition_dimension", &Config::get_partition_dimension, &Config::set_partition_dimension)
      .def_property("subpartition_dimension", &Config::get_subpartition_dimension, &Config::set_subpartition_dimension)
      .def_property("total_iterations", &Config::get_total_iterations, &Config::set_total_iterations)
      .def_property("check_overlapped_walls", &Config::get_check_overlapped_walls, &Config::set_check_overlapped_walls)
      .def_property("reaction_class_cleanup_periodicity", &Config::get_reaction_class_cleanup_periodicity, &Config::set_reaction_class_cleanup_periodicity)
      .def_property("species_cleanup_periodicity", &Config::get_species_cleanup_periodicity, &Config::set_species_cleanup_periodicity)
      .def_property("sort_molecules", &Config::get_sort_molecules, &Config::set_sort_molecules)
      .def_property("memory_limit_gb", &Config::get_memory_limit_gb, &Config::set_memory_limit_gb)
      .def_property("initial_iteration", &Config::get_initial_iteration, &Config::set_initial_iteration)
      .def_property("initial_time", &Config::get_initial_time, &Config::set_initial_time)
      .def_property("rng_state", &Config::get_rng_state, &Config::set_rng_state)
    ;
}

std::string GenConfig::export_to_python(std::ostream& out, PythonExportContext& ctx) const {
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
  if (initial_partition_origin != std::vector<float_t>() && !skip_vectors_export()) {
    ss << ind << "initial_partition_origin = " << export_vec_initial_partition_origin(out, ctx, exported_name) << "," << nl;
  }
  if (partition_dimension != 10) {
    ss << ind << "partition_dimension = " << f_to_str(partition_dimension) << "," << nl;
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
  if (is_set(rng_state)) {
    ss << ind << "rng_state = " << rng_state->export_to_python(out, ctx) << "," << nl;
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

std::string GenConfig::export_vec_initial_partition_origin(std::ostream& out, PythonExportContext& ctx, const std::string& parent_name) const {
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

