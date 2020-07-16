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
#include "gen_config.h"
#include "../api/config.h"

namespace MCell {
namespace API {

void GenConfig::check_semantics() const {
}

bool GenConfig::__eq__(const GenConfig& other) const {
  return
    name == other.name &&
    seed == other.seed &&
    time_step == other.time_step &&
    surface_grid_density == other.surface_grid_density &&
    interaction_radius == other.interaction_radius &&
    vacancy_search_distance == other.vacancy_search_distance &&
    center_molecules_on_grid == other.center_molecules_on_grid &&
    initial_partition_origin == other.initial_partition_origin &&
    partition_dimension == other.partition_dimension &&
    subpartition_dimension == other.subpartition_dimension &&
    total_iterations_hint == other.total_iterations_hint;
}

void GenConfig::set_initialized() {
  initialized = true;
}

void GenConfig::set_all_attributes_as_default_or_unset() {
  class_name = "Config";
  seed = 1;
  time_step = 1e-6;
  surface_grid_density = 10000;
  interaction_radius = FLT_UNSET;
  vacancy_search_distance = 10;
  center_molecules_on_grid = false;
  initial_partition_origin = std::vector<float_t>();
  partition_dimension = 10;
  subpartition_dimension = 0.5;
  total_iterations_hint = 1000000;
}

std::string GenConfig::to_str(const std::string ind) const {
  std::stringstream ss;
  ss << get_object_name() << ": " <<
      "seed=" << seed << ", " <<
      "time_step=" << time_step << ", " <<
      "surface_grid_density=" << surface_grid_density << ", " <<
      "interaction_radius=" << interaction_radius << ", " <<
      "vacancy_search_distance=" << vacancy_search_distance << ", " <<
      "center_molecules_on_grid=" << center_molecules_on_grid << ", " <<
      "initial_partition_origin=" << vec_nonptr_to_str(initial_partition_origin, ind + "  ") << ", " <<
      "partition_dimension=" << partition_dimension << ", " <<
      "subpartition_dimension=" << subpartition_dimension << ", " <<
      "total_iterations_hint=" << total_iterations_hint;
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
            const bool,
            const std::vector<float_t>,
            const float_t,
            const float_t,
            const long
          >(),
          py::arg("seed") = 1,
          py::arg("time_step") = 1e-6,
          py::arg("surface_grid_density") = 10000,
          py::arg("interaction_radius") = FLT_UNSET,
          py::arg("vacancy_search_distance") = 10,
          py::arg("center_molecules_on_grid") = false,
          py::arg("initial_partition_origin") = std::vector<float_t>(),
          py::arg("partition_dimension") = 10,
          py::arg("subpartition_dimension") = 0.5,
          py::arg("total_iterations_hint") = 1000000
      )
      .def("check_semantics", &Config::check_semantics)
      .def("__str__", &Config::to_str, py::arg("ind") = std::string(""))
      .def("dump", &Config::dump)
      .def_property("seed", &Config::get_seed, &Config::set_seed)
      .def_property("time_step", &Config::get_time_step, &Config::set_time_step)
      .def_property("surface_grid_density", &Config::get_surface_grid_density, &Config::set_surface_grid_density)
      .def_property("interaction_radius", &Config::get_interaction_radius, &Config::set_interaction_radius)
      .def_property("vacancy_search_distance", &Config::get_vacancy_search_distance, &Config::set_vacancy_search_distance)
      .def_property("center_molecules_on_grid", &Config::get_center_molecules_on_grid, &Config::set_center_molecules_on_grid)
      .def_property("initial_partition_origin", &Config::get_initial_partition_origin, &Config::set_initial_partition_origin)
      .def_property("partition_dimension", &Config::get_partition_dimension, &Config::set_partition_dimension)
      .def_property("subpartition_dimension", &Config::get_subpartition_dimension, &Config::set_subpartition_dimension)
      .def_property("total_iterations_hint", &Config::get_total_iterations_hint, &Config::set_total_iterations_hint)
    ;
}

} // namespace API
} // namespace MCell

