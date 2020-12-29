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
#include "api/python_export_utils.h"
#include "gen_release_pattern.h"
#include "api/release_pattern.h"

namespace MCell {
namespace API {

void GenReleasePattern::check_semantics() const {
}

void GenReleasePattern::set_initialized() {
  initialized = true;
}

void GenReleasePattern::set_all_attributes_as_default_or_unset() {
  class_name = "ReleasePattern";
  name = STR_UNSET;
  release_interval = TIME_INFINITY;
  train_duration = TIME_INFINITY;
  train_interval = TIME_INFINITY;
  number_of_trains = 1;
}

bool GenReleasePattern::__eq__(const ReleasePattern& other) const {
  return
    name == other.name &&
    release_interval == other.release_interval &&
    train_duration == other.train_duration &&
    train_interval == other.train_interval &&
    number_of_trains == other.number_of_trains;
}

bool GenReleasePattern::eq_nonarray_attributes(const ReleasePattern& other, const bool ignore_name) const {
  return
    (ignore_name || name == other.name) &&
    release_interval == other.release_interval &&
    train_duration == other.train_duration &&
    train_interval == other.train_interval &&
    number_of_trains == other.number_of_trains;
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
      .def("__eq__", &ReleasePattern::__eq__, py::arg("other"))
      .def("dump", &ReleasePattern::dump)
      .def_property("name", &ReleasePattern::get_name, &ReleasePattern::set_name)
      .def_property("release_interval", &ReleasePattern::get_release_interval, &ReleasePattern::set_release_interval)
      .def_property("train_duration", &ReleasePattern::get_train_duration, &ReleasePattern::set_train_duration)
      .def_property("train_interval", &ReleasePattern::get_train_interval, &ReleasePattern::set_train_interval)
      .def_property("number_of_trains", &ReleasePattern::get_number_of_trains, &ReleasePattern::set_number_of_trains)
    ;
}

std::string GenReleasePattern::export_to_python(std::ostream& out, PythonExportContext& ctx) const {
  if (!export_even_if_already_exported() && ctx.already_exported(this)) {
    return ctx.get_exported_name(this);
  }
  std::string exported_name = "release_pattern_" + fix_id(name);
  if (!export_even_if_already_exported()) {
    ctx.add_exported(this, exported_name);

  }
  bool str_export = export_as_string_without_newlines();
  std::string nl = "";
  std::string ind = " ";
  std::stringstream ss;
  if (!str_export) {
    nl = "\n";
    ind = "    ";
    ss << exported_name << " = ";
  }
  ss << "m.ReleasePattern(" << nl;
  if (name != STR_UNSET) {
    ss << ind << "name = " << "'" << name << "'" << "," << nl;
  }
  if (release_interval != TIME_INFINITY) {
    ss << ind << "release_interval = " << release_interval << "," << nl;
  }
  if (train_duration != TIME_INFINITY) {
    ss << ind << "train_duration = " << train_duration << "," << nl;
  }
  if (train_interval != TIME_INFINITY) {
    ss << ind << "train_interval = " << train_interval << "," << nl;
  }
  if (number_of_trains != 1) {
    ss << ind << "number_of_trains = " << number_of_trains << "," << nl;
  }
  ss << ")" << nl << nl;
  if (!str_export) {
    out << ss.str();
    return exported_name;
  }
  else {
    return ss.str();
  }
}

} // namespace API
} // namespace MCell

