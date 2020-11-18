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
#include "gen_release_site.h"
#include "api/release_site.h"
#include "api/complex.h"
#include "api/molecule_release_info.h"
#include "api/region.h"
#include "api/release_pattern.h"

namespace MCell {
namespace API {

void GenReleaseSite::check_semantics() const {
  if (!is_set(name)) {
    throw ValueError("Parameter 'name' must be set.");
  }
}

void GenReleaseSite::set_initialized() {
  if (is_set(complex)) {
    complex->set_initialized();
  }
  vec_set_initialized(molecule_list);
  if (is_set(release_pattern)) {
    release_pattern->set_initialized();
  }
  if (is_set(region)) {
    region->set_initialized();
  }
  initialized = true;
}

void GenReleaseSite::set_all_attributes_as_default_or_unset() {
  class_name = "ReleaseSite";
  name = STR_UNSET;
  complex = nullptr;
  molecule_list = std::vector<std::shared_ptr<MoleculeReleaseInfo>>();
  release_time = 0;
  release_pattern = nullptr;
  shape = Shape::UNSET;
  region = nullptr;
  location = VEC3_UNSET;
  site_diameter = 0;
  site_radius = FLT_UNSET;
  number_to_release = FLT_UNSET;
  density = FLT_UNSET;
  concentration = FLT_UNSET;
  release_probability = FLT_UNSET;
}

bool GenReleaseSite::__eq__(const ReleaseSite& other) const {
  return
    name == other.name &&
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
    vec_ptr_eq(molecule_list, other.molecule_list) &&
    release_time == other.release_time &&
    (
      (is_set(release_pattern)) ?
        (is_set(other.release_pattern) ?
          (release_pattern->__eq__(*other.release_pattern)) : 
          false
        ) :
        (is_set(other.release_pattern) ?
          false :
          true
        )
     )  &&
    shape == other.shape &&
    (
      (is_set(region)) ?
        (is_set(other.region) ?
          (region->__eq__(*other.region)) : 
          false
        ) :
        (is_set(other.region) ?
          false :
          true
        )
     )  &&
    location == other.location &&
    site_diameter == other.site_diameter &&
    site_radius == other.site_radius &&
    number_to_release == other.number_to_release &&
    density == other.density &&
    concentration == other.concentration &&
    release_probability == other.release_probability;
}

bool GenReleaseSite::eq_nonarray_attributes(const ReleaseSite& other) const {
  return
    name == other.name &&
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
    true /*molecule_list*/ &&
    release_time == other.release_time &&
    (
      (is_set(release_pattern)) ?
        (is_set(other.release_pattern) ?
          (release_pattern->__eq__(*other.release_pattern)) : 
          false
        ) :
        (is_set(other.release_pattern) ?
          false :
          true
        )
     )  &&
    shape == other.shape &&
    (
      (is_set(region)) ?
        (is_set(other.region) ?
          (region->__eq__(*other.region)) : 
          false
        ) :
        (is_set(other.region) ?
          false :
          true
        )
     )  &&
    location == other.location &&
    site_diameter == other.site_diameter &&
    site_radius == other.site_radius &&
    number_to_release == other.number_to_release &&
    density == other.density &&
    concentration == other.concentration &&
    release_probability == other.release_probability;
}

std::string GenReleaseSite::to_str(const std::string ind) const {
  std::stringstream ss;
  ss << get_object_name() << ": " <<
      "name=" << name << ", " <<
      "\n" << ind + "  " << "complex=" << "(" << ((complex != nullptr) ? complex->to_str(ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "molecule_list=" << vec_ptr_to_str(molecule_list, ind + "  ") << ", " << "\n" << ind + "  " <<
      "release_time=" << release_time << ", " <<
      "\n" << ind + "  " << "release_pattern=" << "(" << ((release_pattern != nullptr) ? release_pattern->to_str(ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "shape=" << shape << ", " <<
      "\n" << ind + "  " << "region=" << "(" << ((region != nullptr) ? region->to_str(ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "location=" << location << ", " <<
      "site_diameter=" << site_diameter << ", " <<
      "site_radius=" << site_radius << ", " <<
      "number_to_release=" << number_to_release << ", " <<
      "density=" << density << ", " <<
      "concentration=" << concentration << ", " <<
      "release_probability=" << release_probability;
  return ss.str();
}

py::class_<ReleaseSite> define_pybinding_ReleaseSite(py::module& m) {
  return py::class_<ReleaseSite, std::shared_ptr<ReleaseSite>>(m, "ReleaseSite")
      .def(
          py::init<
            const std::string&,
            std::shared_ptr<Complex>,
            const std::vector<std::shared_ptr<MoleculeReleaseInfo>>,
            const float_t,
            std::shared_ptr<ReleasePattern>,
            const Shape,
            std::shared_ptr<Region>,
            const Vec3&,
            const float_t,
            const float_t,
            const float_t,
            const float_t,
            const float_t,
            const float_t
          >(),
          py::arg("name"),
          py::arg("complex") = nullptr,
          py::arg("molecule_list") = std::vector<std::shared_ptr<MoleculeReleaseInfo>>(),
          py::arg("release_time") = 0,
          py::arg("release_pattern") = nullptr,
          py::arg("shape") = Shape::UNSET,
          py::arg("region") = nullptr,
          py::arg("location") = VEC3_UNSET,
          py::arg("site_diameter") = 0,
          py::arg("site_radius") = FLT_UNSET,
          py::arg("number_to_release") = FLT_UNSET,
          py::arg("density") = FLT_UNSET,
          py::arg("concentration") = FLT_UNSET,
          py::arg("release_probability") = FLT_UNSET
      )
      .def("check_semantics", &ReleaseSite::check_semantics)
      .def("__str__", &ReleaseSite::to_str, py::arg("ind") = std::string(""))
      .def("__repr__", &ReleaseSite::to_str, py::arg("ind") = std::string(""))
      .def("__eq__", &ReleaseSite::__eq__, py::arg("other"))
      .def("dump", &ReleaseSite::dump)
      .def_property("name", &ReleaseSite::get_name, &ReleaseSite::set_name)
      .def_property("complex", &ReleaseSite::get_complex, &ReleaseSite::set_complex)
      .def_property("molecule_list", &ReleaseSite::get_molecule_list, &ReleaseSite::set_molecule_list)
      .def_property("release_time", &ReleaseSite::get_release_time, &ReleaseSite::set_release_time)
      .def_property("release_pattern", &ReleaseSite::get_release_pattern, &ReleaseSite::set_release_pattern)
      .def_property("shape", &ReleaseSite::get_shape, &ReleaseSite::set_shape)
      .def_property("region", &ReleaseSite::get_region, &ReleaseSite::set_region)
      .def_property("location", &ReleaseSite::get_location, &ReleaseSite::set_location)
      .def_property("site_diameter", &ReleaseSite::get_site_diameter, &ReleaseSite::set_site_diameter)
      .def_property("site_radius", &ReleaseSite::get_site_radius, &ReleaseSite::set_site_radius)
      .def_property("number_to_release", &ReleaseSite::get_number_to_release, &ReleaseSite::set_number_to_release)
      .def_property("density", &ReleaseSite::get_density, &ReleaseSite::set_density)
      .def_property("concentration", &ReleaseSite::get_concentration, &ReleaseSite::set_concentration)
      .def_property("release_probability", &ReleaseSite::get_release_probability, &ReleaseSite::set_release_probability)
    ;
}

} // namespace API
} // namespace MCell

