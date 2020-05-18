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

#include "../api/common.h"

namespace MCell {
namespace API {

void define_pybinding_constants(py::module& m) {
  m.attr("STATE_UNSET") = py::str("STATE_UNSET");
  m.attr("STATE_UNSET_INT") = py::int_(STATE_UNSET_INT);
  m.attr("BOND_UNBOUND") = py::int_(BOND_UNBOUND);
  m.attr("BOND_BOUND") = py::int_(BOND_BOUND);
  py::enum_<Orientation>(m, "Orientation", py::arithmetic())
    .value("Down", Orientation::Down)
    .value("None", Orientation::None)
    .value("Up", Orientation::Up)
    .value("NotSet", Orientation::NotSet)
    .export_values();
  py::enum_<Notification>(m, "Notification", py::arithmetic())
    .value("None", Notification::None)
    .value("Brief", Notification::Brief)
    .value("Full", Notification::Full)
    .export_values();
  py::enum_<WarningLevel>(m, "WarningLevel", py::arithmetic())
    .value("Ignore", WarningLevel::Ignore)
    .value("Warning", WarningLevel::Warning)
    .value("Error", WarningLevel::Error)
    .export_values();
  py::enum_<VizMode>(m, "VizMode", py::arithmetic())
    .value("Ascii", VizMode::Ascii)
    .value("Cellblender", VizMode::Cellblender)
    .export_values();
  py::enum_<Shape>(m, "Shape", py::arithmetic())
    .value("Unset", Shape::Unset)
    .value("Spherical", Shape::Spherical)
    .export_values();
  py::enum_<ExprNodeType>(m, "ExprNodeType", py::arithmetic())
    .value("Leaf", ExprNodeType::Leaf)
    .value("Add", ExprNodeType::Add)
    .value("Sub", ExprNodeType::Sub)
    .export_values();
}

} // namespace API
} // namespace MCell

