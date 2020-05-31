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
#include "gen_species.h"
#include "../api/species.h"
#include "../api/complex_instance.h"
#include "../api/molecule_instance.h"

namespace MCell {
namespace API {

void GenSpecies::check_semantics() const {
  if (!is_set(name)) {
    throw ValueError("Parameter 'name' must be set.");
  }
}

bool GenSpecies::__eq__(const GenSpecies& other) const {
  return
    name == other.name &&
    name == other.name &&
    diffusion_constant_2d == other.diffusion_constant_2d &&
    diffusion_constant_3d == other.diffusion_constant_3d &&
    vec_ptr_eq(molecule_instances, other.molecule_instances) &&
    orientation == other.orientation;
}

void GenSpecies::set_initialized() {
  vec_set_initialized(molecule_instances);
  initialized = true;
}

std::string GenSpecies::to_str(const std::string ind) const {
  std::stringstream ss;
  ss << get_object_name() << ": " <<
      "name=" << name << ", " <<
      "diffusion_constant_2d=" << diffusion_constant_2d << ", " <<
      "diffusion_constant_3d=" << diffusion_constant_3d << ", " <<
      "\n" << ind + "  " << "molecule_instances=" << vec_ptr_to_str(molecule_instances, ind + "  ") << ", " << "\n" << ind + "  " <<
      "orientation=" << orientation;
  return ss.str();
}

py::class_<Species> define_pybinding_Species(py::module& m) {
  return py::class_<Species, ComplexInstance, std::shared_ptr<Species>>(m, "Species")
      .def(
          py::init<
            const std::string&,
            const float_t,
            const float_t,
            const std::vector<std::shared_ptr<MoleculeInstance>>,
            const Orientation
          >(),
          py::arg("name"),
          py::arg("diffusion_constant_2d") = FLT_UNSET,
          py::arg("diffusion_constant_3d") = FLT_UNSET,
          py::arg("molecule_instances") = std::vector<std::shared_ptr<MoleculeInstance>>(),
          py::arg("orientation") = Orientation::NONE
      )
      .def("check_semantics", &Species::check_semantics)
      .def("__str__", &Species::to_str, py::arg("ind") = std::string(""))
      .def("inst", &Species::inst, py::arg("orientation") = Orientation::NOT_SET)
      .def("dump", &Species::dump)
      .def_property("name", &Species::get_name, &Species::set_name)
      .def_property("diffusion_constant_2d", &Species::get_diffusion_constant_2d, &Species::set_diffusion_constant_2d)
      .def_property("diffusion_constant_3d", &Species::get_diffusion_constant_3d, &Species::set_diffusion_constant_3d)
    ;
}

} // namespace API
} // namespace MCell

