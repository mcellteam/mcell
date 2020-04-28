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
#include "gen_instantiation_data.h"
#include "../api/instantiation_data.h"
#include "../api/geometry_object.h"
#include "../api/release_site.h"
#include "../api/species.h"

namespace MCell {
namespace API {

py::class_<InstantiationData> define_pybinding_InstantiationData(py::module& m) {
  return py::class_<InstantiationData>(m, "InstantiationData")
      .def(
          py::init<
          >()

        )
      .def("add_release_site", &InstantiationData::add_release_site, py::arg("s"))
      .def("find_release_site", &InstantiationData::find_release_site, py::arg("name"))
      .def("add_geometry_object", &InstantiationData::add_geometry_object, py::arg("o"))
      .def("find_geometry_object", &InstantiationData::find_geometry_object, py::arg("name"))
      .def("dump", &InstantiationData::dump)
    ;
}

} // namespace API
} // namespace MCell

