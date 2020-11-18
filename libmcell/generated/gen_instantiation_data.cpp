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
#include "gen_instantiation_data.h"
#include "api/instantiation_data.h"
#include "api/geometry_object.h"
#include "api/region.h"
#include "api/release_site.h"
#include "api/subsystem.h"

namespace MCell {
namespace API {

bool GenInstantiationData::__eq__(const InstantiationData& other) const {
  return
    vec_ptr_eq(release_sites, other.release_sites) &&
    vec_ptr_eq(geometry_objects, other.geometry_objects);
}

bool GenInstantiationData::eq_nonarray_attributes(const InstantiationData& other) const {
  return
    true /*release_sites*/ &&
    true /*geometry_objects*/;
}

std::string GenInstantiationData::to_str(const std::string ind) const {
  std::stringstream ss;
  ss << "InstantiationData" << ": " <<
      "\n" << ind + "  " << "release_sites=" << vec_ptr_to_str(release_sites, ind + "  ") << ", " << "\n" << ind + "  " <<
      "geometry_objects=" << vec_ptr_to_str(geometry_objects, ind + "  ");
  return ss.str();
}

py::class_<InstantiationData> define_pybinding_InstantiationData(py::module& m) {
  return py::class_<InstantiationData, std::shared_ptr<InstantiationData>>(m, "InstantiationData")
      .def(
          py::init<
          >()
      )
      .def("__str__", &InstantiationData::to_str, py::arg("ind") = std::string(""))
      .def("__repr__", &InstantiationData::to_str, py::arg("ind") = std::string(""))
      .def("__eq__", &InstantiationData::__eq__, py::arg("other"))
      .def("add_release_site", &InstantiationData::add_release_site, py::arg("s"))
      .def("find_release_site", &InstantiationData::find_release_site, py::arg("name"))
      .def("add_geometry_object", &InstantiationData::add_geometry_object, py::arg("o"))
      .def("find_geometry_object", &InstantiationData::find_geometry_object, py::arg("name"))
      .def("find_volume_compartment", &InstantiationData::find_volume_compartment, py::arg("name"))
      .def("find_surface_compartment", &InstantiationData::find_surface_compartment, py::arg("name"))
      .def("load_bngl_seed_species", &InstantiationData::load_bngl_seed_species, py::arg("file_name"), py::arg("subsystem"), py::arg("default_release_region") = nullptr, py::arg("parameter_overrides") = std::map<std::string, float_t>())
      .def("dump", &InstantiationData::dump)
      .def_property("release_sites", &InstantiationData::get_release_sites, &InstantiationData::set_release_sites)
      .def_property("geometry_objects", &InstantiationData::get_geometry_objects, &InstantiationData::set_geometry_objects)
    ;
}

} // namespace API
} // namespace MCell

