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

#include "../api/common.h"
#include "../api/species.h"

namespace MCell {
namespace API {

void define_pybinding_constants(py::module& m) {
  m.attr("STATE_UNSET") = py::str(STATE_UNSET);
  m.attr("STATE_UNSET_INT") = py::int_(STATE_UNSET_INT);
  m.attr("BOND_UNBOUND") = py::int_(BOND_UNBOUND);
  m.attr("BOND_BOUND") = py::int_(BOND_BOUND);
  m.attr("PARTITION_EDGE_EXTRA_MARGIN_UM") = py::float_(PARTITION_EDGE_EXTRA_MARGIN_UM);
  m.attr("DEFAULT_COUNT_BUFFER_SIZE") = py::int_(DEFAULT_COUNT_BUFFER_SIZE);
  m.attr("ALL_MOLECULES") = py::str(ALL_MOLECULES);
  m.attr("ALL_VOLUME_MOLECULES") = py::str(ALL_VOLUME_MOLECULES);
  m.attr("ALL_SURFACE_MOLECULES") = py::str(ALL_SURFACE_MOLECULES);
  m.attr("AllMolecules") = py::object(py::cast(AllMolecules));
  m.attr("AllVolumeMolecules") = py::object(py::cast(AllVolumeMolecules));
  m.attr("AllSurfaceMolecules") = py::object(py::cast(AllSurfaceMolecules));
  m.attr("MOLECULE_ID_INVALID") = py::int_(MOLECULE_ID_INVALID);
  m.attr("NUMBER_OF_TRAINS_UNLIMITED") = py::int_(NUMBER_OF_TRAINS_UNLIMITED);
  m.attr("TIME_INFINITY") = py::float_(TIME_INFINITY);
}

void define_pybinding_enums(py::module& m) {
  py::enum_<Orientation>(m, "Orientation", py::arithmetic())
    .value("DOWN", Orientation::DOWN)
    .value("NONE", Orientation::NONE)
    .value("UP", Orientation::UP)
    .value("NOT_SET", Orientation::NOT_SET)
    .value("ANY", Orientation::ANY)
    .export_values();
  py::enum_<Notification>(m, "Notification", py::arithmetic())
    .value("NONE", Notification::NONE)
    .value("BRIEF", Notification::BRIEF)
    .value("FULL", Notification::FULL)
    .export_values();
  py::enum_<WarningLevel>(m, "WarningLevel", py::arithmetic())
    .value("IGNORE", WarningLevel::IGNORE)
    .value("WARNING", WarningLevel::WARNING)
    .value("ERROR", WarningLevel::ERROR)
    .export_values();
  py::enum_<VizMode>(m, "VizMode", py::arithmetic())
    .value("ASCII", VizMode::ASCII)
    .value("CELLBLENDER", VizMode::CELLBLENDER)
    .export_values();
  py::enum_<Shape>(m, "Shape", py::arithmetic())
    .value("UNSET", Shape::UNSET)
    .value("SPHERICAL", Shape::SPHERICAL)
    .value("REGION_EXPR", Shape::REGION_EXPR)
    .export_values();
  py::enum_<SurfacePropertyType>(m, "SurfacePropertyType", py::arithmetic())
    .value("UNSET", SurfacePropertyType::UNSET)
    .value("REFLECTIVE", SurfacePropertyType::REFLECTIVE)
    .value("TRANSPARENT", SurfacePropertyType::TRANSPARENT)
    .value("ABSORPTIVE", SurfacePropertyType::ABSORPTIVE)
    .export_values();
  py::enum_<ExprNodeType>(m, "ExprNodeType", py::arithmetic())
    .value("UNSET", ExprNodeType::UNSET)
    .value("LEAF", ExprNodeType::LEAF)
    .value("ADD", ExprNodeType::ADD)
    .value("SUB", ExprNodeType::SUB)
    .export_values();
  py::enum_<RegionNodeType>(m, "RegionNodeType", py::arithmetic())
    .value("UNSET", RegionNodeType::UNSET)
    .value("LEAF_GEOMETRY_OBJECT", RegionNodeType::LEAF_GEOMETRY_OBJECT)
    .value("LEAF_SURFACE_REGION", RegionNodeType::LEAF_SURFACE_REGION)
    .value("UNION", RegionNodeType::UNION)
    .value("DIFFERENCE", RegionNodeType::DIFFERENCE)
    .value("INTERSECT", RegionNodeType::INTERSECT)
    .export_values();
}

} // namespace API
} // namespace MCell

