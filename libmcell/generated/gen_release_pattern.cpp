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
#include "gen_release_pattern.h"
#include "../api/release_pattern.h"

namespace MCell {
namespace API {

void GenReleasePattern::check_semantics() const {
}

bool GenReleasePattern::__eq__(const GenReleasePattern& other) const {
  return
    name == other.name &&
    name == other.name &&
    release_interval == other.release_interval &&
    train_duration == other.train_duration &&
    train_interval == other.train_interval &&
    number_of_trains == other.number_of_trains;
}

void GenReleasePattern::set_initialized() {
  initialized = true;
}

std::string GenReleasePattern::to_str(const std::string ind) const {
  std::stringstream ss;
  ss << get_object_name() << ": " <<
      "name=" << name << ", " <<
      "release_interval=" << release_interval << ", " <<
      "train_duration=" << train_duration << ", " <<
      "train_interval=" << train_interval << ", " <<
      "number_of_trains=" << number_of_trains;
  return ss.str();
}

py::class_<ReleasePattern> define_pybinding_ReleasePattern(py::module& m) {
  return py::class_<ReleasePattern, std::shared_ptr<ReleasePattern>>(m, "ReleasePattern")
      .def(
          py::init<
            const std::string&,
            const float_t,
            const float_t,
            const float_t,
            const int
          >(),
          py::arg("name") = STR_UNSET,
          py::arg("release_interval") = TIME_INFINITY,
          py::arg("train_duration") = TIME_INFINITY,
          py::arg("train_interval") = TIME_INFINITY,
          py::arg("number_of_trains") = 1
      )
      .def("check_semantics", &ReleasePattern::check_semantics)
      .def("__str__", &ReleasePattern::to_str, py::arg("ind") = std::string(""))
      .def("dump", &ReleasePattern::dump)
      .def_property("name", &ReleasePattern::get_name, &ReleasePattern::set_name)
      .def_property("release_interval", &ReleasePattern::get_release_interval, &ReleasePattern::set_release_interval)
      .def_property("train_duration", &ReleasePattern::get_train_duration, &ReleasePattern::set_train_duration)
      .def_property("train_interval", &ReleasePattern::get_train_interval, &ReleasePattern::set_train_interval)
      .def_property("number_of_trains", &ReleasePattern::get_number_of_trains, &ReleasePattern::set_number_of_trains)
    ;
}

} // namespace API
} // namespace MCell

