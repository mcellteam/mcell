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
#include "gen_complex_instance.h"
#include "../api/complex_instance.h"
#include "../api/molecule_type.h"

namespace MCell {
namespace API {

SemRes GenComplexInstance::check_semantics(std::ostream& out) const {
  return SemRes::OK;
}

std::string GenComplexInstance::to_str() const {
  std::stringstream ss;
  ss << get_object_name() << ": " <<
      "molecule_types=" << vec_ptr_to_str(molecule_types);
  return ss.str();
}

py::class_<ComplexInstance> define_pybinding_ComplexInstance(py::module& m) {
  return py::class_<ComplexInstance, std::shared_ptr<ComplexInstance>>(m, "ComplexInstance")
      .def(
          py::init<
            const std::vector<std::shared_ptr<MoleculeType>>
          >()
,          py::arg("molecule_types") = std::vector<std::shared_ptr<MoleculeType>>()
        )
      .def("check_semantics", &ComplexInstance::check_semantics_cerr)
      .def("__str__", &ComplexInstance::to_str)
      .def("dump", &ComplexInstance::dump)
      .def_property("molecule_types", &ComplexInstance::get_molecule_types, &ComplexInstance::set_molecule_types)
    ;
}

} // namespace API
} // namespace MCell

