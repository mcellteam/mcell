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
#include "gen_model.h"
#include "../api/model.h"
#include "../api/geometry_object.h"
#include "../api/instantiation_data.h"
#include "../api/release_site.h"
#include "../api/species.h"
#include "../api/subsystem.h"

namespace MCell {
namespace API {

py::class_<Model> define_pybinding_Model(py::module& m) {
  return py::class_<Model, std::shared_ptr<Model>>(m, "Model")
      .def(
          py::init<
          >()

        )
      .def("run_iterations", &Model::run_iterations, py::arg("iterations"))
      .def("add_subsystem", &Model::add_subsystem, py::arg("subsystem"))
      .def("add_instantiation_data", &Model::add_instantiation_data, py::arg("instantiation_data"))
      .def("add_species", &Model::add_species, py::arg("s"))
      .def("find_species", &Model::find_species, py::arg("name"))
      .def("add_release_site", &Model::add_release_site, py::arg("s"))
      .def("find_release_site", &Model::find_release_site, py::arg("name"))
      .def("add_geometry_object", &Model::add_geometry_object, py::arg("o"))
      .def("find_geometry_object", &Model::find_geometry_object, py::arg("name"))
      .def("dump", &Model::dump)
    ;
}

} // namespace API
} // namespace MCell

