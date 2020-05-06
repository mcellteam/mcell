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
#include "gen_notifications.h"
#include "../api/notifications.h"

namespace MCell {
namespace API {

void GenNotifications::check_semantics() const {
}

std::string GenNotifications::to_str(const std::string ind) const {
  std::stringstream ss;
  ss << get_object_name() << ": " <<
      "probability_report=" << probability_report << ", " <<
      "diffusion_constant_report=" << diffusion_constant_report << ", " <<
      "final_summary=" << final_summary << ", " <<
      "iteration_report=" << iteration_report << ", " <<
      "varying_probability_report=" << varying_probability_report << ", " <<
      "progress_report=" << progress_report << ", " <<
      "release_event_report=" << release_event_report << ", " <<
      "molecule_collision_report=" << molecule_collision_report;
  return ss.str();
}

py::class_<Notifications> define_pybinding_Notifications(py::module& m) {
  return py::class_<Notifications, std::shared_ptr<Notifications>>(m, "Notifications")
      .def(
          py::init<
            const bool,
            const Notification,
            const bool,
            const bool,
            const bool,
            const bool,
            const bool,
            const bool
          >(),
          py::arg("probability_report") = true,
          py::arg("diffusion_constant_report") = Notification::Brief,
          py::arg("final_summary") = true,
          py::arg("iteration_report") = true,
          py::arg("varying_probability_report") = true,
          py::arg("progress_report") = true,
          py::arg("release_event_report") = true,
          py::arg("molecule_collision_report") = true
      )
      .def("check_semantics", &Notifications::check_semantics)
      .def("__str__", &Notifications::to_str, py::arg("ind") = std::string(""))
      .def("dump", &Notifications::dump)
      .def_property("probability_report", &Notifications::get_probability_report, &Notifications::set_probability_report)
      .def_property("diffusion_constant_report", &Notifications::get_diffusion_constant_report, &Notifications::set_diffusion_constant_report)
      .def_property("final_summary", &Notifications::get_final_summary, &Notifications::set_final_summary)
      .def_property("iteration_report", &Notifications::get_iteration_report, &Notifications::set_iteration_report)
      .def_property("varying_probability_report", &Notifications::get_varying_probability_report, &Notifications::set_varying_probability_report)
      .def_property("progress_report", &Notifications::get_progress_report, &Notifications::set_progress_report)
      .def_property("release_event_report", &Notifications::get_release_event_report, &Notifications::set_release_event_report)
      .def_property("molecule_collision_report", &Notifications::get_molecule_collision_report, &Notifications::set_molecule_collision_report)
    ;
}

} // namespace API
} // namespace MCell

