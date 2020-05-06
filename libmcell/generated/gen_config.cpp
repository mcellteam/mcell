/******************************************************************************
 *
 * Copyright (C) 2020 by
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

#include <sstream>
#include <pybind11/stl.h>
#include "gen_config.h"
#include "../api/config.h"

namespace MCell {
namespace API {

void GenConfig::check_semantics() const {
}

std::string GenConfig::to_str(const std::string ind) const {
  std::stringstream ss;
  ss << get_object_name() << ": " <<
      "time_step=" << time_step << ", " <<
      "surface_grid_density=" << surface_grid_density << ", " <<
      "accurate_3d_reactions=" << accurate_3d_reactions << ", " <<
      "center_molecules_on_grid=" << center_molecules_on_grid << ", " <<
      "microscopic_reversibility=" << microscopic_reversibility;
  return ss.str();
}

py::class_<Config> define_pybinding_Config(py::module& m) {
  return py::class_<Config, std::shared_ptr<Config>>(m, "Config")
      .def(
          py::init<
            const float_t,
            const float_t,
            const bool,
            const bool,
            const bool
          >(),
          py::arg("time_step") = 1e-6,
          py::arg("surface_grid_density") = 10000,
          py::arg("accurate_3d_reactions") = true,
          py::arg("center_molecules_on_grid") = false,
          py::arg("microscopic_reversibility") = false
      )
      .def("check_semantics", &Config::check_semantics)
      .def("__str__", &Config::to_str, py::arg("ind") = std::string(""))
      .def("dump", &Config::dump)
      .def_property("time_step", &Config::get_time_step, &Config::set_time_step)
      .def_property("surface_grid_density", &Config::get_surface_grid_density, &Config::set_surface_grid_density)
      .def_property("accurate_3d_reactions", &Config::get_accurate_3d_reactions, &Config::set_accurate_3d_reactions)
      .def_property("center_molecules_on_grid", &Config::get_center_molecules_on_grid, &Config::set_center_molecules_on_grid)
      .def_property("microscopic_reversibility", &Config::get_microscopic_reversibility, &Config::set_microscopic_reversibility)
    ;
}

} // namespace API
} // namespace MCell

