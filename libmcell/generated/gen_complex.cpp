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
#include "gen_complex.h"
#include "api/complex.h"
#include "api/elementary_molecule_instance.h"
#include "api/species.h"

namespace MCell {
namespace API {

void GenComplex::check_semantics() const {
}

void GenComplex::set_initialized() {
  vec_set_initialized(elementary_molecule_instances);
  initialized = true;
}

void GenComplex::set_all_attributes_as_default_or_unset() {
  class_name = "Complex";
  name = STR_UNSET;
  elementary_molecule_instances = std::vector<std::shared_ptr<ElementaryMoleculeInstance>>();
  orientation = Orientation::DEFAULT;
  compartment_name = STR_UNSET;
}

bool GenComplex::__eq__(const Complex& other) const {
  return
    name == other.name &&
    name == other.name &&
    vec_ptr_eq(elementary_molecule_instances, other.elementary_molecule_instances) &&
    orientation == other.orientation &&
    compartment_name == other.compartment_name;
}

std::string GenComplex::to_str(const std::string ind) const {
  std::stringstream ss;
  ss << get_object_name() << ": " <<
      "name=" << name << ", " <<
      "\n" << ind + "  " << "elementary_molecule_instances=" << vec_ptr_to_str(elementary_molecule_instances, ind + "  ") << ", " << "\n" << ind + "  " <<
      "orientation=" << orientation << ", " <<
      "compartment_name=" << compartment_name;
  return ss.str();
}

py::class_<Complex> define_pybinding_Complex(py::module& m) {
  return py::class_<Complex, std::shared_ptr<Complex>>(m, "Complex")
      .def(
          py::init<
            const std::string&,
            const std::vector<std::shared_ptr<ElementaryMoleculeInstance>>,
            const Orientation,
            const std::string&
          >(),
          py::arg("name") = STR_UNSET,
          py::arg("elementary_molecule_instances") = std::vector<std::shared_ptr<ElementaryMoleculeInstance>>(),
          py::arg("orientation") = Orientation::DEFAULT,
          py::arg("compartment_name") = STR_UNSET
      )
      .def("check_semantics", &Complex::check_semantics)
      .def("__str__", &Complex::to_str, py::arg("ind") = std::string(""))
      .def("__repr__", &Complex::to_str, py::arg("ind") = std::string(""))
      .def("__eq__", &Complex::__eq__, py::arg("other"))
      .def("to_bngl_str", &Complex::to_bngl_str)
      .def("as_species", &Complex::as_species)
      .def("dump", &Complex::dump)
      .def_property("name", &Complex::get_name, &Complex::set_name)
      .def_property("elementary_molecule_instances", &Complex::get_elementary_molecule_instances, &Complex::set_elementary_molecule_instances)
      .def_property("orientation", &Complex::get_orientation, &Complex::set_orientation)
      .def_property("compartment_name", &Complex::get_compartment_name, &Complex::set_compartment_name)
    ;
}

} // namespace API
} // namespace MCell

