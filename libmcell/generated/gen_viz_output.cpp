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
#include <pybind11/stl.h>
#include "gen_viz_output.h"
#include "../api/viz_output.h"
#include "../api/species.h"

namespace MCell {
namespace API {

void GenVizOutput::check_semantics() const {
  if (!is_set(filename_prefix)) {
    throw ValueError("Parameter 'filename_prefix' must be set.");
  }
}

bool GenVizOutput::__eq__(const GenVizOutput& other) const {
  return
    name == other.name &&
    filename_prefix == other.filename_prefix &&
    vec_ptr_eq(species_list, other.species_list) &&
    mode == other.mode &&
    every_n_timesteps == other.every_n_timesteps;
}

void GenVizOutput::set_initialized() {
  vec_set_initialized(species_list);
  initialized = true;
}

std::string GenVizOutput::to_str(const std::string ind) const {
  std::stringstream ss;
  ss << get_object_name() << ": " <<
      "filename_prefix=" << filename_prefix << ", " <<
      "\n" << ind + "  " << "species_list=" << vec_ptr_to_str(species_list, ind + "  ") << ", " << "\n" << ind + "  " <<
      "mode=" << mode << ", " <<
      "every_n_timesteps=" << every_n_timesteps;
  return ss.str();
}

py::class_<VizOutput> define_pybinding_VizOutput(py::module& m) {
  return py::class_<VizOutput, std::shared_ptr<VizOutput>>(m, "VizOutput")
      .def(
          py::init<
            const std::string&,
            const std::vector<std::shared_ptr<Species>>,
            const VizMode,
            const int
          >(),
          py::arg("filename_prefix"),
          py::arg("species_list") = std::vector<std::shared_ptr<Species>>(),
          py::arg("mode") = VizMode::Ascii,
          py::arg("every_n_timesteps") = 1
      )
      .def("check_semantics", &VizOutput::check_semantics)
      .def("__str__", &VizOutput::to_str, py::arg("ind") = std::string(""))
      .def("dump", &VizOutput::dump)
      .def_property("filename_prefix", &VizOutput::get_filename_prefix, &VizOutput::set_filename_prefix)
      .def_property("species_list", &VizOutput::get_species_list, &VizOutput::set_species_list)
      .def_property("mode", &VizOutput::get_mode, &VizOutput::set_mode)
      .def_property("every_n_timesteps", &VizOutput::get_every_n_timesteps, &VizOutput::set_every_n_timesteps)
    ;
}

} // namespace API
} // namespace MCell

