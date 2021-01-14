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
#include "api/pybind11_stl_include.h"
#include "gen_molecule_release_info.h"
#include "api/molecule_release_info.h"
#include "api/complex.h"

namespace MCell {
namespace API {

void GenMoleculeReleaseInfo::check_semantics() const {
  if (!is_set(complex)) {
    throw ValueError("Parameter 'complex' must be set.");
  }
  if (!is_set(location)) {
    throw ValueError("Parameter 'location' must be set.");
  }
}

void GenMoleculeReleaseInfo::set_initialized() {
  if (is_set(complex)) {
    complex->set_initialized();
  }
  initialized = true;
}

void GenMoleculeReleaseInfo::set_all_attributes_as_default_or_unset() {
  class_name = "MoleculeReleaseInfo";
  complex = nullptr;
  location = std::vector<float_t>();
}

bool GenMoleculeReleaseInfo::__eq__(const MoleculeReleaseInfo& other) const {
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
    location == other.location;
}

bool GenMoleculeReleaseInfo::eq_nonarray_attributes(const MoleculeReleaseInfo& other, const bool ignore_name) const {
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
    true /*location*/;
}

std::string GenMoleculeReleaseInfo::to_str(const std::string ind) const {
  std::stringstream ss;
  ss << get_object_name() << ": " <<
      "\n" << ind + "  " << "complex=" << "(" << ((complex != nullptr) ? complex->to_str(ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "location=" << vec_nonptr_to_str(location, ind + "  ");
  return ss.str();
}

py::class_<MoleculeReleaseInfo> define_pybinding_MoleculeReleaseInfo(py::module& m) {
  return py::class_<MoleculeReleaseInfo, std::shared_ptr<MoleculeReleaseInfo>>(m, "MoleculeReleaseInfo")
      .def(
          py::init<
            std::shared_ptr<Complex>,
            const std::vector<float_t>
          >(),
          py::arg("complex"),
          py::arg("location")
      )
      .def("check_semantics", &MoleculeReleaseInfo::check_semantics)
      .def("__str__", &MoleculeReleaseInfo::to_str, py::arg("ind") = std::string(""))
      .def("__eq__", &MoleculeReleaseInfo::__eq__, py::arg("other"))
      .def("dump", &MoleculeReleaseInfo::dump)
      .def_property("complex", &MoleculeReleaseInfo::get_complex, &MoleculeReleaseInfo::set_complex)
      .def_property("location", &MoleculeReleaseInfo::get_location, &MoleculeReleaseInfo::set_location)
    ;
}

} // namespace API
} // namespace MCell

