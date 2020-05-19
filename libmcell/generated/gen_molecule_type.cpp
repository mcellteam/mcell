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
#include "gen_molecule_type.h"
#include "../api/molecule_type.h"
#include "../api/component_instance.h"
#include "../api/component_type.h"
#include "../api/molecule_instance.h"

namespace MCell {
namespace API {

void GenMoleculeType::check_semantics() const {
  if (!is_set(name)) {
    throw ValueError("Parameter 'name' must be set.");
  }
}

std::string GenMoleculeType::to_str(const std::string ind) const {
  std::stringstream ss;
  ss << get_object_name() << ": " <<
      "name=" << name << ", " <<
      "\n" << ind + "  " << "components=" << vec_ptr_to_str(components, ind + "  ") << ", " << "\n" << ind + "  " <<
      "diffusion_constant_2d=" << diffusion_constant_2d << ", " <<
      "diffusion_constant_3d=" << diffusion_constant_3d;
  return ss.str();
}

bool GenMoleculeType::__eq__(const GenMoleculeType& other) const {
  return
    name == other.name &&
    name == other.name &&
    vec_ptr_eq(components, other.components) &&
    diffusion_constant_2d == other.diffusion_constant_2d &&
    diffusion_constant_3d == other.diffusion_constant_3d;
}

py::class_<MoleculeType> define_pybinding_MoleculeType(py::module& m) {
  return py::class_<MoleculeType, std::shared_ptr<MoleculeType>>(m, "MoleculeType")
      .def(
          py::init<
            const std::string&,
            const std::vector<std::shared_ptr<ComponentType>>,
            const float_t,
            const float_t
          >(),
          py::arg("name"),
          py::arg("components") = std::vector<std::shared_ptr<ComponentType>>(),
          py::arg("diffusion_constant_2d") = FLT_UNSET,
          py::arg("diffusion_constant_3d") = FLT_UNSET
      )
      .def("check_semantics", &MoleculeType::check_semantics)
      .def("__str__", &MoleculeType::to_str, py::arg("ind") = std::string(""))
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

