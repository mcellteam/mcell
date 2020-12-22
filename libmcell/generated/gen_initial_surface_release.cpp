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
#include "gen_initial_surface_release.h"
#include "api/initial_surface_release.h"
#include "api/complex.h"

namespace MCell {
namespace API {

void GenInitialSurfaceRelease::check_semantics() const {
  if (!is_set(complex)) {
    throw ValueError("Parameter 'complex' must be set.");
  }
}

void GenInitialSurfaceRelease::set_initialized() {
  if (is_set(complex)) {
    complex->set_initialized();
  }
  initialized = true;
}

void GenInitialSurfaceRelease::set_all_attributes_as_default_or_unset() {
  class_name = "InitialSurfaceRelease";
  complex = nullptr;
  number_to_release = INT_UNSET;
  density = FLT_UNSET;
}

bool GenInitialSurfaceRelease::__eq__(const InitialSurfaceRelease& other) const {
  return
    (
      (is_set(complex)) ?
        (is_set(other.complex) ?
          (complex->__eq__(*other.complex)) : 
          false
        ) :
        (is_set(other.complex) ?
          false :
          true
        )
     )  &&
    number_to_release == other.number_to_release &&
    density == other.density;
}

bool GenInitialSurfaceRelease::eq_nonarray_attributes(const InitialSurfaceRelease& other, const bool ignore_name) const {
  return
    (
      (is_set(complex)) ?
        (is_set(other.complex) ?
          (complex->__eq__(*other.complex)) : 
          false
        ) :
        (is_set(other.complex) ?
          false :
          true
        )
     )  &&
    number_to_release == other.number_to_release &&
    density == other.density;
}

std::string GenInitialSurfaceRelease::to_str(const std::string ind) const {
  std::stringstream ss;
  ss << get_object_name() << ": " <<
      "\n" << ind + "  " << "complex=" << "(" << ((complex != nullptr) ? complex->to_str(ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "number_to_release=" << number_to_release << ", " <<
      "density=" << density;
  return ss.str();
}

py::class_<InitialSurfaceRelease> define_pybinding_InitialSurfaceRelease(py::module& m) {
  return py::class_<InitialSurfaceRelease, std::shared_ptr<InitialSurfaceRelease>>(m, "InitialSurfaceRelease")
      .def(
          py::init<
            std::shared_ptr<Complex>,
            const int,
            const float_t
          >(),
          py::arg("complex"),
          py::arg("number_to_release") = INT_UNSET,
          py::arg("density") = FLT_UNSET
      )
      .def("check_semantics", &InitialSurfaceRelease::check_semantics)
      .def("__str__", &InitialSurfaceRelease::to_str, py::arg("ind") = std::string(""))
      .def("__eq__", &InitialSurfaceRelease::__eq__, py::arg("other"))
      .def("dump", &InitialSurfaceRelease::dump)
      .def_property("complex", &InitialSurfaceRelease::get_complex, &InitialSurfaceRelease::set_complex)
      .def_property("number_to_release", &InitialSurfaceRelease::get_number_to_release, &InitialSurfaceRelease::set_number_to_release)
      .def_property("density", &InitialSurfaceRelease::get_density, &InitialSurfaceRelease::set_density)
    ;
}

std::string GenInitialSurfaceRelease::export_to_python(std::ostream& out) const {
  std::string name = "TODO";
  std::stringstream ss;
  ss << name << " = GenInitialSurfaceRelease(\n";
  ss << "  complex = " << complex->export_to_python(out) << ",\n";
  if (number_to_release != INT_UNSET) {
    ss << "  number_to_release = " << number_to_release << ",\n";
  }
  if (density != FLT_UNSET) {
    ss << "  density = " << density << ",\n";
  }
  ss << ")\n\n";
  out << ss.str();
  return name;
}

} // namespace API
} // namespace MCell

