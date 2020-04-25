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
#include "gen_partition.h"
#include "../api/partition.h"

namespace MCell {
namespace API {

SemRes GenPartition::check_semantics(std::ostream& out) const {
  return SemRes::OK;
}

std::string GenPartition::to_str(const std::string ind) const {
  std::stringstream ss;
  ss << get_object_name() << ": " <<
      "parition_dimension=" << parition_dimension << ", " <<
      "subparition_dimension=" << subparition_dimension;
  return ss.str();
}

py::class_<Partition> define_pybinding_Partition(py::module& m) {
  return py::class_<Partition, std::shared_ptr<Partition>>(m, "Partition")
      .def(
          py::init<
            const float_t,
            const float_t
          >(),
          py::arg("parition_dimension") = 10,
          py::arg("subparition_dimension") = 0.5
        )
      .def("check_semantics", &Partition::check_semantics_cerr)
      .def("__str__", &Partition::to_str, py::arg("ind") = std::string(""))
      .def("dump", &Partition::dump)
      .def_property("parition_dimension", &Partition::get_parition_dimension, &Partition::set_parition_dimension)
      .def_property("subparition_dimension", &Partition::get_subparition_dimension, &Partition::set_subparition_dimension)
    ;
}

} // namespace API
} // namespace MCell

