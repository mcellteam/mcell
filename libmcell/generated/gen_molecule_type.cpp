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
#include "gen_molecule_type.h"
#include "../api/molecule_type.h"
#include "../api/component_instance.h"
#include "../api/component_type.h"
#include "../api/molecule_instance.h"

namespace MCell {
namespace API {

SemRes GenMoleculeType::check_semantics(std::ostream& out) const {
  if (!is_set(name)) {
    out << get_object_name() << ": Parameter 'name' must be set.\n";
    return SemRes::ERROR;
  }
  return SemRes::OK;
}

std::string GenMoleculeType::to_str() const {
  std::stringstream ss;
  ss << get_object_name() << ": " <<
      "name=" << name << ", " <<
      "components=" << vec_ptr_to_str(components) << ", " <<
      "diffusion_constant_2d=" << diffusion_constant_2d << ", " <<
      "diffusion_constant_3d=" << diffusion_constant_3d;
  return ss.str();
}

py::class_<MoleculeType> define_pybinding_MoleculeType(py::module& m) {
  return py::class_<MoleculeType>(m, "MoleculeType")
      .def(
          py::init<
            const std::string&,
            const std::vector<ComponentType*>,
            const float_t,
            const float_t
          >()
,          py::arg("name"),
          py::arg("components") = std::vector<ComponentType*>(),
          py::arg("diffusion_constant_2d") = FLT_UNSET,
          py::arg("diffusion_constant_3d") = FLT_UNSET
        )
      .def("check_semantics", &MoleculeType::check_semantics_cerr)
      .def("__str__", &MoleculeType::to_str)
      .def("inst", &MoleculeType::inst, py::arg("components"))
      .def("dump", &MoleculeType::dump)
      .def_property("name", &MoleculeType::get_name, &MoleculeType::set_name)
      .def_property("components", &MoleculeType::get_components, &MoleculeType::set_components)
      .def_property("diffusion_constant_2d", &MoleculeType::get_diffusion_constant_2d, &MoleculeType::set_diffusion_constant_2d)
      .def_property("diffusion_constant_3d", &MoleculeType::get_diffusion_constant_3d, &MoleculeType::set_diffusion_constant_3d)
    ;
}

} // namespace API
} // namespace MCell

