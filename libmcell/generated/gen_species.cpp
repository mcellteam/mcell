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
#include "gen_species.h"
#include "../api/species.h"
#include "../api/complex.h"
#include "../api/elementary_molecule_instance.h"
#include "../api/species.h"

namespace MCell {
namespace API {

void GenSpecies::check_semantics() const {
}

bool GenSpecies::__eq__(const GenSpecies& other) const {
  return
    name == other.name &&
    name == other.name &&
    diffusion_constant_2d == other.diffusion_constant_2d &&
    diffusion_constant_3d == other.diffusion_constant_3d &&
    custom_time_step == other.custom_time_step &&
    custom_space_step == other.custom_space_step &&
    target_only == other.target_only &&
    name == other.name &&
    vec_ptr_eq(elementary_molecule_instances, other.elementary_molecule_instances) &&
    orientation == other.orientation &&
    compartment_name == other.compartment_name;
}

void GenSpecies::set_initialized() {
  vec_set_initialized(elementary_molecule_instances);
  initialized = true;
}

void GenSpecies::set_all_attributes_as_default_or_unset() {
  class_name = "Species";
  name = STR_UNSET;
  diffusion_constant_2d = FLT_UNSET;
  diffusion_constant_3d = FLT_UNSET;
  custom_time_step = FLT_UNSET;
  custom_space_step = FLT_UNSET;
  target_only = false;
  name = STR_UNSET;
  elementary_molecule_instances = std::vector<std::shared_ptr<ElementaryMoleculeInstance>>();
  orientation = Orientation::DEFAULT;
  compartment_name = STR_UNSET;
}

std::string GenSpecies::to_str(const std::string ind) const {
  std::stringstream ss;
  ss << get_object_name() << ": " <<
      "name=" << name << ", " <<
      "diffusion_constant_2d=" << diffusion_constant_2d << ", " <<
      "diffusion_constant_3d=" << diffusion_constant_3d << ", " <<
      "custom_time_step=" << custom_time_step << ", " <<
      "custom_space_step=" << custom_space_step << ", " <<
      "target_only=" << target_only << ", " <<
      "name=" << name << ", " <<
      "\n" << ind + "  " << "elementary_molecule_instances=" << vec_ptr_to_str(elementary_molecule_instances, ind + "  ") << ", " << "\n" << ind + "  " <<
      "orientation=" << orientation << ", " <<
      "compartment_name=" << compartment_name;
  return ss.str();
}

py::class_<Species> define_pybinding_Species(py::module& m) {
  return py::class_<Species, Complex, std::shared_ptr<Species>>(m, "Species")
      .def(
          py::init<
            const std::string&,
            const float_t,
            const float_t,
            const float_t,
            const float_t,
            const bool,
            const std::vector<std::shared_ptr<ElementaryMoleculeInstance>>,
            const Orientation,
            const std::string&
          >(),
          py::arg("name") = STR_UNSET,
          py::arg("diffusion_constant_2d") = FLT_UNSET,
          py::arg("diffusion_constant_3d") = FLT_UNSET,
          py::arg("custom_time_step") = FLT_UNSET,
          py::arg("custom_space_step") = FLT_UNSET,
          py::arg("target_only") = false,
          py::arg("elementary_molecule_instances") = std::vector<std::shared_ptr<ElementaryMoleculeInstance>>(),
          py::arg("orientation") = Orientation::DEFAULT,
          py::arg("compartment_name") = STR_UNSET
      )
      .def("check_semantics", &Species::check_semantics)
      .def("__str__", &Species::to_str, py::arg("ind") = std::string(""))
      .def("inst", &Species::inst, py::arg("orientation") = Orientation::DEFAULT, py::arg("compartment_name") = STR_UNSET)
      .def("dump", &Species::dump)
      .def_property("name", &Species::get_name, &Species::set_name)
      .def_property("diffusion_constant_2d", &Species::get_diffusion_constant_2d, &Species::set_diffusion_constant_2d)
      .def_property("diffusion_constant_3d", &Species::get_diffusion_constant_3d, &Species::set_diffusion_constant_3d)
      .def_property("custom_time_step", &Species::get_custom_time_step, &Species::set_custom_time_step)
      .def_property("custom_space_step", &Species::get_custom_space_step, &Species::set_custom_space_step)
      .def_property("target_only", &Species::get_target_only, &Species::set_target_only)
    ;
}

} // namespace API
} // namespace MCell

