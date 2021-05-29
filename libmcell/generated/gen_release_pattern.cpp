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
#include "api/pybind11_stl_include.h"
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

std::shared_ptr<ReleasePattern> GenReleasePattern::copy_release_pattern() const {
  if (initialized) {
    throw RuntimeError("Object of class ReleasePattern cannot be cloned with 'copy' after this object was used in model initialization.");
  }

  std::shared_ptr<ReleasePattern> res = std::make_shared<ReleasePattern>(DefaultCtorArgType());
  res->class_name = class_name;
  res->name = name;
  res->release_interval = release_interval;
  res->train_duration = train_duration;
  res->train_interval = train_interval;
  res->number_of_trains = number_of_trains;

  return res;
}

std::shared_ptr<ReleasePattern> GenReleasePattern::deepcopy_release_pattern(py::dict) const {
  if (initialized) {
    throw RuntimeError("Object of class ReleasePattern cannot be cloned with 'deepcopy' after this object was used in model initialization.");
  }

  std::shared_ptr<ReleasePattern> res = std::make_shared<ReleasePattern>(DefaultCtorArgType());
  res->class_name = class_name;
  res->name = name;
  res->release_interval = release_interval;
  res->train_duration = train_duration;
  res->train_interval = train_interval;
  res->number_of_trains = number_of_trains;

  return res;
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

std::string GenReleasePattern::to_str(const bool all_details, const std::string ind) const {
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
  return py::class_<ReleasePattern, std::shared_ptr<ReleasePattern>>(m, "ReleasePattern", "Defines a release pattern that specifies repeating molecule releases. \nCan be used by a ReleaseSite.\n")
      .def(
          py::init<
            const std::string&,
            const double,
            const double,
            const double,
            const int
          >(),
          py::arg("name") = STR_UNSET,
          py::arg("release_interval") = TIME_INFINITY,
          py::arg("train_duration") = TIME_INFINITY,
          py::arg("train_interval") = TIME_INFINITY,
          py::arg("number_of_trains") = 1
      )
      .def("check_semantics", &ReleasePattern::check_semantics)
      .def("__copy__", &ReleasePattern::copy_release_pattern)
      .def("__deepcopy__", &ReleasePattern::deepcopy_release_pattern, py::arg("memo"))
      .def("__str__", &ReleasePattern::to_str, py::arg("all_details") = false, py::arg("ind") = std::string(""))
      .def("__eq__", &ReleasePattern::__eq__, py::arg("other"))
      .def("dump", &ReleasePattern::dump)
      .def_property("name", &ReleasePattern::get_name, &ReleasePattern::set_name, "Name of the release pattern.")
      .def_property("release_interval", &ReleasePattern::get_release_interval, &ReleasePattern::set_release_interval, "During a train of releases, release molecules after every t seconds. \nDefault is to release only once.\n")
      .def_property("train_duration", &ReleasePattern::get_train_duration, &ReleasePattern::set_train_duration, "The train of releases lasts for t seconds before turning off. \nDefault is to never turn off.\n")
      .def_property("train_interval", &ReleasePattern::get_train_interval, &ReleasePattern::set_train_interval, "A new train of releases happens every t seconds. \nDefault is to never have a new train. \nThe train interval must not be shorter than the train duration.\n")
      .def_property("number_of_trains", &ReleasePattern::get_number_of_trains, &ReleasePattern::set_number_of_trains, "Repeat the release process for n trains of releases. Default is one train.\nFor unlimited number of trains use a constant NUMBER_OF_TRAINS_UNLIMITED.\n")
    ;
}

std::string GenReleasePattern::export_to_python(std::ostream& out, PythonExportContext& ctx) {
  if (!export_even_if_already_exported() && ctx.already_exported(this)) {
    return ctx.get_exported_name(this);
  }
  std::string exported_name = std::string("release_pattern") + "_" + (is_set(name) ? fix_id(name) : std::to_string(ctx.postinc_counter("release_pattern")));
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
    ss << ind << "release_interval = " << f_to_str(release_interval) << "," << nl;
  }
  if (train_duration != TIME_INFINITY) {
    ss << ind << "train_duration = " << f_to_str(train_duration) << "," << nl;
  }
  if (train_interval != TIME_INFINITY) {
    ss << ind << "train_interval = " << f_to_str(train_interval) << "," << nl;
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

