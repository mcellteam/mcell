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
#include "../api/mcell.h"
namespace MCell {
namespace API {

py::class_<Model> define_pybinding_Model(py::module& m) {
  return py::class_<Model>(m, "Model")
      .def(
          py::init<
          >()

        )
      .def("run_iterations", &Model::run_iterations, py::arg("iterations"))
      .def("add_subsystem", &Model::add_subsystem, py::arg("subsystem"))
      .def("add_instantiation_data", &Model::add_instantiation_data, py::arg("instantiation_data"))
    ;
}

} // namespace API
} // namespace MCell

