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
  high_reaction_probability = WarningLevel::IGNORE;
}

bool GenWarnings::__eq__(const Warnings& other) const {
  return
    high_reaction_probability == other.high_reaction_probability;
}

bool GenWarnings::eq_nonarray_attributes(const Warnings& other, const bool ignore_name) const {
  return
    high_reaction_probability == other.high_reaction_probability;
}

std::string GenWarnings::to_str(const std::string ind) const {
  std::stringstream ss;
  ss << get_object_name() << ": " <<
      "high_reaction_probability=" << high_reaction_probability;
  return ss.str();
}

py::class_<Warnings> define_pybinding_Warnings(py::module& m) {
  return py::class_<Warnings, std::shared_ptr<Warnings>>(m, "Warnings", "This class contains warnings settings. For now it contains only one configurable \nwarning.\n")
      .def(
          py::init<
            const WarningLevel
          >(),
          py::arg("high_reaction_probability") = WarningLevel::IGNORE
      )
      .def("check_semantics", &Warnings::check_semantics)
      .def("__str__", &Warnings::to_str, py::arg("ind") = std::string(""))
      .def("__eq__", &Warnings::__eq__, py::arg("other"))
      .def("dump", &Warnings::dump)
      .def_property("high_reaction_probability", &Warnings::get_high_reaction_probability, &Warnings::set_high_reaction_probability, "Print a warning when a bimolecular reaction probability is over 0.5 but less or equal than 1.\nWarning when probability is greater than 1 is always printed.\nCannot be set to WarningLevel.ERROR.\n")
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
  if (high_reaction_probability != WarningLevel::IGNORE) {
    ss << ind << "high_reaction_probability = " << high_reaction_probability << "," << nl;
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

