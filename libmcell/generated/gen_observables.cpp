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
#include "gen_observables.h"
#include "../api/observables.h"
#include "../api/count.h"
#include "../api/viz_output.h"

namespace MCell {
namespace API {

std::string GenObservables::to_str(const std::string ind) const {
  std::stringstream ss;
  ss << "Observables" << ": " <<
      "\n" << ind + "  " << "viz_outputs=" << vec_ptr_to_str(viz_outputs, ind + "  ") << ", " << "\n" << ind + "  " <<
      "counts=" << vec_ptr_to_str(counts, ind + "  ");
  return ss.str();
}

py::class_<Observables> define_pybinding_Observables(py::module& m) {
  return py::class_<Observables, std::shared_ptr<Observables>>(m, "Observables")
      .def(
          py::init<
          >()
      )
      .def("__str__", &Observables::to_str, py::arg("ind") = std::string(""))
      .def("add_viz_output", &Observables::add_viz_output, py::arg("viz_output"))
      .def("add_count", &Observables::add_count, py::arg("count"))
      .def("dump", &Observables::dump)
      .def_property("viz_outputs", &Observables::get_viz_outputs, &Observables::set_viz_outputs)
      .def_property("counts", &Observables::get_counts, &Observables::set_counts)
    ;
}

} // namespace API
} // namespace MCell
