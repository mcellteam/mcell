/******************************************************************************
 *
 * Copyright (C) 2021 by
 * The Salk Institute for Biological Studies
 *
 * Use of this source code is governed by an MIT-style
 * license that can be found in the LICENSE file or at
 * https://opensource.org/licenses/MIT.
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
  molecule_placement_failure = WarningLevel::ERROR;
}

std::shared_ptr<Warnings> GenWarnings::copy_warnings() const {
  std::shared_ptr<Warnings> res = std::make_shared<Warnings>(DefaultCtorArgType());
  res->class_name = class_name;
  res->high_reaction_probability = high_reaction_probability;
  res->molecule_placement_failure = molecule_placement_failure;

  return res;
}

std::shared_ptr<Warnings> GenWarnings::deepcopy_warnings(py::dict) const {
  std::shared_ptr<Warnings> res = std::make_shared<Warnings>(DefaultCtorArgType());
  res->class_name = class_name;
  res->high_reaction_probability = high_reaction_probability;
  res->molecule_placement_failure = molecule_placement_failure;

  return res;
}

bool GenWarnings::__eq__(const Warnings& other) const {
  return
    high_reaction_probability == other.high_reaction_probability &&
    molecule_placement_failure == other.molecule_placement_failure;
}

bool GenWarnings::eq_nonarray_attributes(const Warnings& other, const bool ignore_name) const {
  return
    high_reaction_probability == other.high_reaction_probability &&
    molecule_placement_failure == other.molecule_placement_failure;
}

std::string GenWarnings::to_str(const bool all_details, const std::string ind) const {
  std::stringstream ss;
  ss << get_object_name() << ": " <<
      "high_reaction_probability=" << high_reaction_probability << ", " <<
      "molecule_placement_failure=" << molecule_placement_failure;
  return ss.str();
}

py::class_<Warnings> define_pybinding_Warnings(py::module& m) {
  return py::class_<Warnings, std::shared_ptr<Warnings>>(m, "Warnings", "This class contains warnings settings. For now it contains only one configurable \nwarning.\n")
      .def(
          py::init<
            const WarningLevel,
            const WarningLevel
          >(),
          py::arg("high_reaction_probability") = WarningLevel::IGNORE,
          py::arg("molecule_placement_failure") = WarningLevel::ERROR
      )
      .def("check_semantics", &Warnings::check_semantics)
      .def("__copy__", &Warnings::copy_warnings)
      .def("__deepcopy__", &Warnings::deepcopy_warnings, py::arg("memo"))
      .def("__str__", &Warnings::to_str, py::arg("all_details") = false, py::arg("ind") = std::string(""))
      .def("__eq__", &Warnings::__eq__, py::arg("other"))
      .def("dump", &Warnings::dump)
      .def_property("high_reaction_probability", &Warnings::get_high_reaction_probability, &Warnings::set_high_reaction_probability, "Print a warning when a bimolecular reaction probability is over 0.5 but less or equal than 1.\nWarning when probability is greater than 1 is always printed.\nCannot be set to WarningLevel.ERROR.\n")
      .def_property("molecule_placement_failure", &Warnings::get_molecule_placement_failure, &Warnings::set_molecule_placement_failure, "Print a warning or end with an error when a release of a molecule fails.\n")
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
  if (molecule_placement_failure != WarningLevel::ERROR) {
    ss << ind << "molecule_placement_failure = " << molecule_placement_failure << "," << nl;
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

