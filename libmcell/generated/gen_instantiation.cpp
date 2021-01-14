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
#include "gen_instantiation.h"
#include "api/instantiation.h"
#include "api/geometry_object.h"
#include "api/region.h"
#include "api/release_site.h"
#include "api/subsystem.h"

namespace MCell {
namespace API {

bool GenInstantiation::__eq__(const Instantiation& other) const {
  return
    vec_ptr_eq(release_sites, other.release_sites) &&
    vec_ptr_eq(geometry_objects, other.geometry_objects);
}

bool GenInstantiation::eq_nonarray_attributes(const Instantiation& other, const bool ignore_name) const {
  return
    true /*release_sites*/ &&
    true /*geometry_objects*/;
}

std::string GenInstantiation::to_str(const std::string ind) const {
  std::stringstream ss;
  ss << "Instantiation" << ": " <<
      "\n" << ind + "  " << "release_sites=" << vec_ptr_to_str(release_sites, ind + "  ") << ", " << "\n" << ind + "  " <<
      "geometry_objects=" << vec_ptr_to_str(geometry_objects, ind + "  ");
  return ss.str();
}

py::class_<Instantiation> define_pybinding_Instantiation(py::module& m) {
  return py::class_<Instantiation, std::shared_ptr<Instantiation>>(m, "Instantiation")
      .def(
          py::init<
          >()
      )
      .def("__str__", &Instantiation::to_str, py::arg("ind") = std::string(""))
      .def("__eq__", &Instantiation::__eq__, py::arg("other"))
      .def("add_release_site", &Instantiation::add_release_site, py::arg("s"))
      .def("find_release_site", &Instantiation::find_release_site, py::arg("name"))
      .def("add_geometry_object", &Instantiation::add_geometry_object, py::arg("o"))
      .def("find_geometry_object", &Instantiation::find_geometry_object, py::arg("name"))
      .def("find_volume_compartment", &Instantiation::find_volume_compartment, py::arg("name"))
      .def("find_surface_compartment", &Instantiation::find_surface_compartment, py::arg("name"))
      .def("load_bngl_seed_species", &Instantiation::load_bngl_seed_species, py::arg("file_name"), py::arg("subsystem"), py::arg("default_release_region") = nullptr, py::arg("parameter_overrides") = std::map<std::string, float_t>())
      .def("dump", &Instantiation::dump)
      .def_property("release_sites", &Instantiation::get_release_sites, &Instantiation::set_release_sites)
      .def_property("geometry_objects", &Instantiation::get_geometry_objects, &Instantiation::set_geometry_objects)
    ;
}

} // namespace API
} // namespace MCell

