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
#include "api/python_export.h"
#include "gen_observables.h"
#include "api/observables.h"
#include "api/count.h"
#include "api/subsystem.h"
#include "api/viz_output.h"

namespace MCell {
namespace API {

bool GenObservables::__eq__(const Observables& other) const {
  return
    vec_ptr_eq(viz_outputs, other.viz_outputs) &&
    vec_ptr_eq(counts, other.counts);
}

bool GenObservables::eq_nonarray_attributes(const Observables& other, const bool ignore_name) const {
  return
    true /*viz_outputs*/ &&
    true /*counts*/;
}

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
      .def("__eq__", &Observables::__eq__, py::arg("other"))
      .def("add_viz_output", &Observables::add_viz_output, py::arg("viz_output"))
      .def("add_count", &Observables::add_count, py::arg("count"))
      .def("find_count", &Observables::find_count, py::arg("name"))
      .def("load_bngl_observables", &Observables::load_bngl_observables, py::arg("file_name"), py::arg("subsystem"), py::arg("output_files_prefix") = "", py::arg("parameter_overrides") = std::map<std::string, float_t>())
      .def("dump", &Observables::dump)
      .def_property("viz_outputs", &Observables::get_viz_outputs, &Observables::set_viz_outputs)
      .def_property("counts", &Observables::get_counts, &Observables::set_counts)
    ;
}

} // namespace API
} // namespace MCell

