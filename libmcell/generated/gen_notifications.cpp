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
#include "gen_notifications.h"
#include "api/notifications.h"

namespace MCell {
namespace API {

void GenNotifications::check_semantics() const {
}

void GenNotifications::set_initialized() {
  initialized = true;
}

void GenNotifications::set_all_attributes_as_default_or_unset() {
  class_name = "Notifications";
  bng_verbosity_level = 0;
  rxn_and_species_report = false;
  simulation_stats_every_n_iterations = 0;
  rxn_probability_changed = true;
}

std::shared_ptr<Notifications> GenNotifications::copy_notifications() const {
  if (initialized) {
    throw RuntimeError("Object of class Notifications cannot be cloned with 'copy' after this object was used in model initialization.");
  }

  std::shared_ptr<Notifications> res = std::make_shared<Notifications>(DefaultCtorArgType());
  res->class_name = class_name;
  res->bng_verbosity_level = bng_verbosity_level;
  res->rxn_and_species_report = rxn_and_species_report;
  res->simulation_stats_every_n_iterations = simulation_stats_every_n_iterations;
  res->rxn_probability_changed = rxn_probability_changed;

  return res;
}

std::shared_ptr<Notifications> GenNotifications::deepcopy_notifications(py::dict) const {
  if (initialized) {
    throw RuntimeError("Object of class Notifications cannot be cloned with 'deepcopy' after this object was used in model initialization.");
  }

  std::shared_ptr<Notifications> res = std::make_shared<Notifications>(DefaultCtorArgType());
  res->class_name = class_name;
  res->bng_verbosity_level = bng_verbosity_level;
  res->rxn_and_species_report = rxn_and_species_report;
  res->simulation_stats_every_n_iterations = simulation_stats_every_n_iterations;
  res->rxn_probability_changed = rxn_probability_changed;

  return res;
}

bool GenNotifications::__eq__(const Notifications& other) const {
  return
    bng_verbosity_level == other.bng_verbosity_level &&
    rxn_and_species_report == other.rxn_and_species_report &&
    simulation_stats_every_n_iterations == other.simulation_stats_every_n_iterations &&
    rxn_probability_changed == other.rxn_probability_changed;
}

bool GenNotifications::eq_nonarray_attributes(const Notifications& other, const bool ignore_name) const {
  return
    bng_verbosity_level == other.bng_verbosity_level &&
    rxn_and_species_report == other.rxn_and_species_report &&
    simulation_stats_every_n_iterations == other.simulation_stats_every_n_iterations &&
    rxn_probability_changed == other.rxn_probability_changed;
}

std::string GenNotifications::to_str(const std::string ind) const {
  std::stringstream ss;
  ss << get_object_name() << ": " <<
      "bng_verbosity_level=" << bng_verbosity_level << ", " <<
      "rxn_and_species_report=" << rxn_and_species_report << ", " <<
      "simulation_stats_every_n_iterations=" << simulation_stats_every_n_iterations << ", " <<
      "rxn_probability_changed=" << rxn_probability_changed;
  return ss.str();
}

py::class_<Notifications> define_pybinding_Notifications(py::module& m) {
  return py::class_<Notifications, std::shared_ptr<Notifications>>(m, "Notifications")
      .def(
          py::init<
            const int,
            const bool,
            const int,
            const bool
          >(),
          py::arg("bng_verbosity_level") = 0,
          py::arg("rxn_and_species_report") = false,
          py::arg("simulation_stats_every_n_iterations") = 0,
          py::arg("rxn_probability_changed") = true
      )
      .def("check_semantics", &Notifications::check_semantics)
      .def("__copy__", &Notifications::copy_notifications)
      .def("__deepcopy__", &Notifications::deepcopy_notifications, py::arg("memo"))
      .def("__str__", &Notifications::to_str, py::arg("ind") = std::string(""))
      .def("__eq__", &Notifications::__eq__, py::arg("other"))
      .def("dump", &Notifications::dump)
      .def_property("bng_verbosity_level", &Notifications::get_bng_verbosity_level, &Notifications::set_bng_verbosity_level, "Sets verbosity level that enables printouts of extra information on BioNetGen \nspecies and rules created and used during simulation.\n")
      .def_property("rxn_and_species_report", &Notifications::get_rxn_and_species_report, &Notifications::set_rxn_and_species_report, "When set to True, simulation generates files rxn_report_SEED.txt, and \nspecies_report_SEED.txt that contain details on reaction classes and species \nthat were created based on reaction rules.   \n")
      .def_property("simulation_stats_every_n_iterations", &Notifications::get_simulation_stats_every_n_iterations, &Notifications::set_simulation_stats_every_n_iterations, "When set to a value other than 0, internal simulation stats will be printed. \n")
      .def_property("rxn_probability_changed", &Notifications::get_rxn_probability_changed, &Notifications::set_rxn_probability_changed, "When True, information that a reaction's probability has changed is printed during simulation.    \n")
    ;
}

std::string GenNotifications::export_to_python(std::ostream& out, PythonExportContext& ctx) {
  if (!export_even_if_already_exported() && ctx.already_exported(this)) {
    return ctx.get_exported_name(this);
  }
  std::string exported_name = "notifications_" + std::to_string(ctx.postinc_counter("notifications"));
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
  ss << "m.Notifications(" << nl;
  if (bng_verbosity_level != 0) {
    ss << ind << "bng_verbosity_level = " << bng_verbosity_level << "," << nl;
  }
  if (rxn_and_species_report != false) {
    ss << ind << "rxn_and_species_report = " << rxn_and_species_report << "," << nl;
  }
  if (simulation_stats_every_n_iterations != 0) {
    ss << ind << "simulation_stats_every_n_iterations = " << simulation_stats_every_n_iterations << "," << nl;
  }
  if (rxn_probability_changed != true) {
    ss << ind << "rxn_probability_changed = " << rxn_probability_changed << "," << nl;
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

