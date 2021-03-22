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
  return py::class_<Warnings, std::shared_ptr<Warnings>>(m, "Warnings", "This is a placeholder for future warnings settings. Empty for now.")
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

std::string GenWarnings::export_to_python(std::ostream& out, PythonExportContext& ctx) {
  if (!export_even_if_already_exported() && ctx.already_exported(this)) {
    return ctx.get_exported_name(this);
  }
  std::string exported_name = "warnings_" + std::to_string(ctx.postinc_counter("warnings"));
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
  ss << "m.Warnings(" << nl;
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

