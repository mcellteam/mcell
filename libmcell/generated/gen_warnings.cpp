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
#include "../api/warnings.h"

namespace MCell {
namespace API {

void GenWarnings::check_semantics() const {
}

bool GenWarnings::__eq__(const GenWarnings& other) const {
  return
    name == other.name &&
    molecule_collision_report == other.molecule_collision_report &&
    degenerate_polygons == other.degenerate_polygons &&
    negative_diffusion_constant == other.negative_diffusion_constant &&
    missing_surface_orientation == other.missing_surface_orientation &&
    negative_reaction_rate == other.negative_reaction_rate &&
    useless_volume_orientation == other.useless_volume_orientation &&
    high_reaction_probability == other.high_reaction_probability &&
    lifetime_too_short == other.lifetime_too_short &&
    lifetime_threshold == other.lifetime_threshold &&
    missed_reactions == other.missed_reactions &&
    missed_reactions_threshold == other.missed_reactions_threshold;
}

void GenWarnings::set_initialized() {
  initialized = true;
}

void GenWarnings::set_all_attributes_as_default_or_unset() {
  class_name = "Warnings";
  molecule_collision_report = WarningLevel::WARNING;
  degenerate_polygons = WarningLevel::WARNING;
  negative_diffusion_constant = WarningLevel::WARNING;
  missing_surface_orientation = WarningLevel::ERROR;
  negative_reaction_rate = WarningLevel::WARNING;
  useless_volume_orientation = WarningLevel::WARNING;
  high_reaction_probability = WarningLevel::IGNORE;
  lifetime_too_short = WarningLevel::WARNING;
  lifetime_threshold = 50;
  missed_reactions = WarningLevel::WARNING;
  missed_reactions_threshold = 0.00100000004749745;
}

std::string GenWarnings::to_str(const std::string ind) const {
  std::stringstream ss;
  ss << get_object_name() << ": " <<
      "molecule_collision_report=" << molecule_collision_report << ", " <<
      "degenerate_polygons=" << degenerate_polygons << ", " <<
      "negative_diffusion_constant=" << negative_diffusion_constant << ", " <<
      "missing_surface_orientation=" << missing_surface_orientation << ", " <<
      "negative_reaction_rate=" << negative_reaction_rate << ", " <<
      "useless_volume_orientation=" << useless_volume_orientation << ", " <<
      "high_reaction_probability=" << high_reaction_probability << ", " <<
      "lifetime_too_short=" << lifetime_too_short << ", " <<
      "lifetime_threshold=" << lifetime_threshold << ", " <<
      "missed_reactions=" << missed_reactions << ", " <<
      "missed_reactions_threshold=" << missed_reactions_threshold;
  return ss.str();
}

py::class_<Warnings> define_pybinding_Warnings(py::module& m) {
  return py::class_<Warnings, std::shared_ptr<Warnings>>(m, "Warnings")
      .def(
          py::init<
            const WarningLevel,
            const WarningLevel,
            const WarningLevel,
            const WarningLevel,
            const WarningLevel,
            const WarningLevel,
            const WarningLevel,
            const WarningLevel,
            const float_t,
            const WarningLevel,
            const float_t
          >(),
          py::arg("molecule_collision_report") = WarningLevel::WARNING,
          py::arg("degenerate_polygons") = WarningLevel::WARNING,
          py::arg("negative_diffusion_constant") = WarningLevel::WARNING,
          py::arg("missing_surface_orientation") = WarningLevel::ERROR,
          py::arg("negative_reaction_rate") = WarningLevel::WARNING,
          py::arg("useless_volume_orientation") = WarningLevel::WARNING,
          py::arg("high_reaction_probability") = WarningLevel::IGNORE,
          py::arg("lifetime_too_short") = WarningLevel::WARNING,
          py::arg("lifetime_threshold") = 50,
          py::arg("missed_reactions") = WarningLevel::WARNING,
          py::arg("missed_reactions_threshold") = 0.00100000004749745
      )
      .def("check_semantics", &Warnings::check_semantics)
      .def("__str__", &Warnings::to_str, py::arg("ind") = std::string(""))
      .def("dump", &Warnings::dump)
      .def_property("molecule_collision_report", &Warnings::get_molecule_collision_report, &Warnings::set_molecule_collision_report)
      .def_property("degenerate_polygons", &Warnings::get_degenerate_polygons, &Warnings::set_degenerate_polygons)
      .def_property("negative_diffusion_constant", &Warnings::get_negative_diffusion_constant, &Warnings::set_negative_diffusion_constant)
      .def_property("missing_surface_orientation", &Warnings::get_missing_surface_orientation, &Warnings::set_missing_surface_orientation)
      .def_property("negative_reaction_rate", &Warnings::get_negative_reaction_rate, &Warnings::set_negative_reaction_rate)
      .def_property("useless_volume_orientation", &Warnings::get_useless_volume_orientation, &Warnings::set_useless_volume_orientation)
      .def_property("high_reaction_probability", &Warnings::get_high_reaction_probability, &Warnings::set_high_reaction_probability)
      .def_property("lifetime_too_short", &Warnings::get_lifetime_too_short, &Warnings::set_lifetime_too_short)
      .def_property("lifetime_threshold", &Warnings::get_lifetime_threshold, &Warnings::set_lifetime_threshold)
      .def_property("missed_reactions", &Warnings::get_missed_reactions, &Warnings::set_missed_reactions)
      .def_property("missed_reactions_threshold", &Warnings::get_missed_reactions_threshold, &Warnings::set_missed_reactions_threshold)
    ;
}

} // namespace API
} // namespace MCell

