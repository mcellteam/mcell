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
#include "gen_release_site.h"
#include "../api/release_site.h"
#include "../api/region.h"
#include "../api/species.h"

namespace MCell {
namespace API {

void GenReleaseSite::check_semantics() const {
  if (!is_set(name)) {
    throw ValueError("Parameter 'name' must be set.");
  }
  if (!is_set(species)) {
    throw ValueError("Parameter 'species' must be set.");
  }
}

bool GenReleaseSite::__eq__(const GenReleaseSite& other) const {
  return
    name == other.name &&
    name == other.name &&
    (
      (species != nullptr) ?
        ( (other.species != nullptr) ?
          (species->__eq__(*other.species)) : 
          false
        ) :
        ( (other.species != nullptr) ?
          false :
          true
        )
     )  &&
    orientation == other.orientation &&
    shape == other.shape &&
    (
      (region != nullptr) ?
        ( (other.region != nullptr) ?
          (region->__eq__(*other.region)) : 
          false
        ) :
        ( (other.region != nullptr) ?
          false :
          true
        )
     )  &&
    location == other.location &&
    site_diameter == other.site_diameter &&
    site_radius == other.site_radius &&
    number_to_release == other.number_to_release &&
    density == other.density &&
    release_probability == other.release_probability;
}

void GenReleaseSite::set_initialized() {
  if (is_set(species)) {
    species->set_initialized();
  }
  if (is_set(region)) {
    region->set_initialized();
  }
  initialized = true;
}

std::string GenReleaseSite::to_str(const std::string ind) const {
  std::stringstream ss;
  ss << get_object_name() << ": " <<
      "name=" << name << ", " <<
      "\n" << ind + "  " << "species=" << "(" << ((species != nullptr) ? species->to_str(ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "orientation=" << orientation << ", " <<
      "shape=" << shape << ", " <<
      "\n" << ind + "  " << "region=" << "(" << ((region != nullptr) ? region->to_str(ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "location=" << location << ", " <<
      "site_diameter=" << site_diameter << ", " <<
      "site_radius=" << site_radius << ", " <<
      "number_to_release=" << number_to_release << ", " <<
      "density=" << density << ", " <<
      "release_probability=" << release_probability;
  return ss.str();
}

py::class_<ReleaseSite> define_pybinding_ReleaseSite(py::module& m) {
  return py::class_<ReleaseSite, std::shared_ptr<ReleaseSite>>(m, "ReleaseSite")
      .def(
          py::init<
            const std::string&,
            std::shared_ptr<Species>,
            const Orientation,
            const Shape,
            std::shared_ptr<Region>,
            const Vec3&,
            const float_t,
            const float_t,
            const int,
            const float_t,
            const float_t
          >(),
          py::arg("name"),
          py::arg("species"),
          py::arg("orientation") = Orientation::NONE,
          py::arg("shape") = Shape::UNSET,
          py::arg("region") = nullptr,
          py::arg("location") = VEC3_UNSET,
          py::arg("site_diameter") = 0,
          py::arg("site_radius") = FLT_UNSET,
          py::arg("number_to_release") = INT_UNSET,
          py::arg("density") = FLT_UNSET,
          py::arg("release_probability") = FLT_UNSET
      )
      .def("check_semantics", &ReleaseSite::check_semantics)
      .def("__str__", &ReleaseSite::to_str, py::arg("ind") = std::string(""))
      .def("dump", &ReleaseSite::dump)
      .def_property("name", &ReleaseSite::get_name, &ReleaseSite::set_name)
      .def_property("species", &ReleaseSite::get_species, &ReleaseSite::set_species)
      .def_property("orientation", &ReleaseSite::get_orientation, &ReleaseSite::set_orientation)
      .def_property("shape", &ReleaseSite::get_shape, &ReleaseSite::set_shape)
      .def_property("region", &ReleaseSite::get_region, &ReleaseSite::set_region)
      .def_property("location", &ReleaseSite::get_location, &ReleaseSite::set_location)
      .def_property("site_diameter", &ReleaseSite::get_site_diameter, &ReleaseSite::set_site_diameter)
      .def_property("site_radius", &ReleaseSite::get_site_radius, &ReleaseSite::set_site_radius)
      .def_property("number_to_release", &ReleaseSite::get_number_to_release, &ReleaseSite::set_number_to_release)
      .def_property("density", &ReleaseSite::get_density, &ReleaseSite::set_density)
      .def_property("release_probability", &ReleaseSite::get_release_probability, &ReleaseSite::set_release_probability)
    ;
}

} // namespace API
} // namespace MCell

