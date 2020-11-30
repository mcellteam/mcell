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
#include "gen_warnings.h"
#include "api/warnings.h"

namespace MCell {
namespace API {

void GenWarnings::check_semantics() const {
}

void GenWarnings::set_initialized() {
  initialized = true;
}

void GenWarnings::set_all_attributes_as_default_or_unset() {
  class_name = "Warnings";
}

bool GenWarnings::__eq__(const Warnings& other) const {
  return
true ;
}

bool GenWarnings::eq_nonarray_attributes(const Warnings& other, const bool ignore_name) const {
  return
true ;
}

std::string GenWarnings::to_str(const std::string ind) const {
  std::stringstream ss;
  ss << get_object_name();
  return ss.str();
}

py::class_<Warnings> define_pybinding_Warnings(py::module& m) {
  return py::class_<Warnings, std::shared_ptr<Warnings>>(m, "Warnings")
      .def(
          py::init<
          >()

      )
      .def("check_semantics", &Warnings::check_semantics)
      .def("__str__", &Warnings::to_str, py::arg("ind") = std::string(""))
      .def("__eq__", &Warnings::__eq__, py::arg("other"))
      .def("dump", &Warnings::dump)
    ;
}

} // namespace API
} // namespace MCell

