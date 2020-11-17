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
#include "gen_viz_output.h"
#include "api/viz_output.h"
#include "api/species.h"

namespace MCell {
namespace API {

void GenVizOutput::check_semantics() const {
  if (!is_set(output_files_prefix)) {
    throw ValueError("Parameter 'output_files_prefix' must be set.");
  }
}

void GenVizOutput::set_initialized() {
  vec_set_initialized(species_list);
  initialized = true;
}

void GenVizOutput::set_all_attributes_as_default_or_unset() {
  class_name = "VizOutput";
  output_files_prefix = STR_UNSET;
  species_list = std::vector<std::shared_ptr<Species>>();
  all_species = false;
  mode = VizMode::ASCII;
  every_n_timesteps = 1;
}

bool GenVizOutput::__eq__(const VizOutput& other) const {
  return
    name == other.name &&
    output_files_prefix == other.output_files_prefix &&
    vec_ptr_eq(species_list, other.species_list) &&
    all_species == other.all_species &&
    mode == other.mode &&
    every_n_timesteps == other.every_n_timesteps;
}

std::string GenVizOutput::to_str(const std::string ind) const {
  std::stringstream ss;
  ss << get_object_name() << ": " <<
      "output_files_prefix=" << output_files_prefix << ", " <<
      "\n" << ind + "  " << "species_list=" << vec_ptr_to_str(species_list, ind + "  ") << ", " << "\n" << ind + "  " <<
      "all_species=" << all_species << ", " <<
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
            const bool,
            const VizMode,
            const float_t
          >(),
          py::arg("output_files_prefix"),
          py::arg("species_list") = std::vector<std::shared_ptr<Species>>(),
          py::arg("all_species") = false,
          py::arg("mode") = VizMode::ASCII,
          py::arg("every_n_timesteps") = 1
      )
      .def("check_semantics", &VizOutput::check_semantics)
      .def("__str__", &VizOutput::to_str, py::arg("ind") = std::string(""))
      .def("__repr__", &VizOutput::to_str, py::arg("ind") = std::string(""))
      .def("__eq__", &VizOutput::__eq__, py::arg("other"))
      .def("dump", &VizOutput::dump)
      .def_property("output_files_prefix", &VizOutput::get_output_files_prefix, &VizOutput::set_output_files_prefix)
      .def_property("species_list", &VizOutput::get_species_list, &VizOutput::set_species_list)
      .def_property("all_species", &VizOutput::get_all_species, &VizOutput::set_all_species)
      .def_property("mode", &VizOutput::get_mode, &VizOutput::set_mode)
      .def_property("every_n_timesteps", &VizOutput::get_every_n_timesteps, &VizOutput::set_every_n_timesteps)
    ;
}

} // namespace API
} // namespace MCell

