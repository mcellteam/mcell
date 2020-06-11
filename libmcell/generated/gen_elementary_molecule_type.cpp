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
#include "gen_elementary_molecule_type.h"
#include "../api/elementary_molecule_type.h"
#include "../api/component_instance.h"
#include "../api/component_type.h"
#include "../api/elementary_molecule_instance.h"

namespace MCell {
namespace API {

void GenElementaryMoleculeType::check_semantics() const {
  if (!is_set(name)) {
    throw ValueError("Parameter 'name' must be set.");
  }
}

bool GenElementaryMoleculeType::__eq__(const GenElementaryMoleculeType& other) const {
  return
    name == other.name &&
    name == other.name &&
    vec_ptr_eq(components, other.components) &&
    diffusion_constant_2d == other.diffusion_constant_2d &&
    diffusion_constant_3d == other.diffusion_constant_3d;
}

void GenElementaryMoleculeType::set_initialized() {
  vec_set_initialized(components);
  initialized = true;
}

std::string GenElementaryMoleculeType::to_str(const std::string ind) const {
  std::stringstream ss;
  ss << get_object_name() << ": " <<
      "name=" << name << ", " <<
      "\n" << ind + "  " << "components=" << vec_ptr_to_str(components, ind + "  ") << ", " << "\n" << ind + "  " <<
      "diffusion_constant_2d=" << diffusion_constant_2d << ", " <<
      "diffusion_constant_3d=" << diffusion_constant_3d;
  return ss.str();
}

py::class_<ElementaryMoleculeType> define_pybinding_ElementaryMoleculeType(py::module& m) {
  return py::class_<ElementaryMoleculeType, std::shared_ptr<ElementaryMoleculeType>>(m, "ElementaryMoleculeType")
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
      .def("check_semantics", &ElementaryMoleculeType::check_semantics)
      .def("__str__", &ElementaryMoleculeType::to_str, py::arg("ind") = std::string(""))
      .def("inst", &ElementaryMoleculeType::inst, py::arg("components"))
      .def("dump", &ElementaryMoleculeType::dump)
      .def_property("name", &ElementaryMoleculeType::get_name, &ElementaryMoleculeType::set_name)
      .def_property("components", &ElementaryMoleculeType::get_components, &ElementaryMoleculeType::set_components)
      .def_property("diffusion_constant_2d", &ElementaryMoleculeType::get_diffusion_constant_2d, &ElementaryMoleculeType::set_diffusion_constant_2d)
      .def_property("diffusion_constant_3d", &ElementaryMoleculeType::get_diffusion_constant_3d, &ElementaryMoleculeType::set_diffusion_constant_3d)
    ;
}

} // namespace API
} // namespace MCell

