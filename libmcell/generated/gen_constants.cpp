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

#include "api/api_common.h"
#include "api/species.h"

namespace MCell {
namespace API {

void define_pybinding_constants(py::module& m) {
  m.attr("STATE_UNSET") = py::str(STATE_UNSET);
  m.attr("STATE_UNSET_INT") = py::int_(STATE_UNSET_INT);
  m.attr("BOND_UNBOUND") = py::int_(BOND_UNBOUND);
  m.attr("BOND_BOUND") = py::int_(BOND_BOUND);
  m.attr("BOND_ANY") = py::int_(BOND_ANY);
  m.attr("PARTITION_EDGE_EXTRA_MARGIN_UM") = py::float_(PARTITION_EDGE_EXTRA_MARGIN_UM);
  m.attr("DEFAULT_COUNT_BUFFER_SIZE") = py::int_(DEFAULT_COUNT_BUFFER_SIZE);
  m.attr("ALL_MOLECULES") = py::str(ALL_MOLECULES);
  m.attr("ALL_VOLUME_MOLECULES") = py::str(ALL_VOLUME_MOLECULES);
  m.attr("ALL_SURFACE_MOLECULES") = py::str(ALL_SURFACE_MOLECULES);
  m.attr("DEFAULT_CHECKPOINTS_DIR") = py::str(DEFAULT_CHECKPOINTS_DIR);
  m.attr("DEFAULT_SEED_DIR_PREFIX") = py::str(DEFAULT_SEED_DIR_PREFIX);
  m.attr("DEFAULT_SEED_DIR_DIGITS") = py::int_(DEFAULT_SEED_DIR_DIGITS);
  m.attr("DEFAULT_ITERATION_DIR_PREFIX") = py::str(DEFAULT_ITERATION_DIR_PREFIX);
  m.attr("AllMolecules") = py::object(py::cast(AllMolecules));
  m.attr("AllVolumeMolecules") = py::object(py::cast(AllVolumeMolecules));
  m.attr("AllSurfaceMolecules") = py::object(py::cast(AllSurfaceMolecules));
  m.attr("ID_INVALID") = py::int_(ID_INVALID);
  m.attr("NUMBER_OF_TRAINS_UNLIMITED") = py::int_(NUMBER_OF_TRAINS_UNLIMITED);
  m.attr("TIME_INFINITY") = py::float_(TIME_INFINITY);
  m.attr("INT_UNSET") = py::int_(INT_UNSET);
  m.attr("FLT_UNSET") = py::float_(FLT_UNSET);
  m.attr("RNG_SIZE") = py::int_(RNG_SIZE);
}

void define_pybinding_enums(py::module& m) {
  py::enum_<Orientation>(m, "Orientation", py::arithmetic(), "Orientation of a Complex.\n- DOWN\n\n- NONE\n\n- UP\n\n- NOT_SET\n\n- ANY\n\n- DEFAULT: Value DEFAULT means NONE for volume complexes and UP for surface complexes.\n")
    .value("DOWN", Orientation::DOWN)
    .value("NONE", Orientation::NONE)
    .value("UP", Orientation::UP)
    .value("NOT_SET", Orientation::NOT_SET)
    .value("ANY", Orientation::ANY)
    .value("DEFAULT", Orientation::DEFAULT)
    .export_values();
  py::enum_<Notification>(m, "Notification", py::arithmetic(), "- NONE\n\n- BRIEF\n\n- FULL\n\n")
    .value("NONE", Notification::NONE)
    .value("BRIEF", Notification::BRIEF)
    .value("FULL", Notification::FULL)
    .export_values();
  py::enum_<WarningLevel>(m, "WarningLevel", py::arithmetic(), "- IGNORE: Do something sensible and continue silently.\n- WARNING: Do something sensible but emit a warning message.\n- ERROR: Treat the warning as an error and stop.\n")
    .value("IGNORE", WarningLevel::IGNORE)
    .value("WARNING", WarningLevel::WARNING)
    .value("ERROR", WarningLevel::ERROR)
    .export_values();
  py::enum_<VizMode>(m, "VizMode", py::arithmetic(), "- ASCII: Readable molecule visualization output.\n- CELLBLENDER_V1: Binary molecule visualization output used by MCell3, format v1.\nAllows only limited length of species name (256 chars) and \ndoes not contain molecule IDs.   \n\n- CELLBLENDER: Binary molecule visualization output, format v2.\n")
    .value("ASCII", VizMode::ASCII)
    .value("CELLBLENDER_V1", VizMode::CELLBLENDER_V1)
    .value("CELLBLENDER", VizMode::CELLBLENDER)
    .export_values();
  py::enum_<Shape>(m, "Shape", py::arithmetic(), "- UNSET\n\n- SPHERICAL\n\n- REGION_EXPR\n\n- LIST\n\n- COMPARTMENT\n\n")
    .value("UNSET", Shape::UNSET)
    .value("SPHERICAL", Shape::SPHERICAL)
    .value("REGION_EXPR", Shape::REGION_EXPR)
    .value("LIST", Shape::LIST)
    .value("COMPARTMENT", Shape::COMPARTMENT)
    .export_values();
  py::enum_<SurfacePropertyType>(m, "SurfacePropertyType", py::arithmetic(), "- UNSET\n\n- REACTIVE: This surface class does not do anything by itself, but it can be used as a reactant in \nreaction rules. \n\n- REFLECTIVE: If used as a surface property for a volume molecule it is reflected by any surface with\nthis surface class. This is the default behavior for volume molecules.\nIf used for a surface molecule it is reflected by the border of the\nsurface with this surface class. \nSetting orientation in affected_complex_pattern allows selective reflection of volume \nmolecules from only the front or back of a surface or selective reflection of surface \nmolecules with only a certain orientation from the surface’s border. \nUsing m.ALL_MOLECULES as affected_complex_pattern has the effect that all \nvolume molecules are reflected by surfaces with this surface class and all surface molecules \nare reflected by the border of the surfaces with this surface class. \nUsing m.ALL_VOLUME_MOLECULES as affected_complex_pattern has the effect that all\nvolume molecules are reflected by surfaces with this surface class. \nUsing m.ALL_SURFACE_MOLECULES as affected_complex_pattern has the effect that all\nsurface molecules are reflected by the border of the surface with this surface class.\n\n- TRANSPARENT: If used as a surface property for a volume molecule it passes through all surfaces with\nthis surface class.  \nIf used for a surface molecule it passes through the border of the surface with this surface \nclass. This is the default behavior for surface molecules.\nSetting orientation in affected_complex_pattern allows selective transparency of volume \nmolecules from only the front or back of a surface or selective transparency for surface \nmolecules with only a certain orientation from the surface’s border. \nTo make a surface with this surface class transparent to all volume molecules,\nuse m.ALL_VOLUME_MOLECULES for affected_complex_pattern. \nTo make a border of the surface with this surface class transparent to all surface molecules,\nuse m.ALL_SURFACE_MOLECULES for the affected_complex_pattern. \nUsing m.ALL_MOLECULES for affected_complex_pattern has the effect that surfaces with this surface class \nare transparent to all volume molecules and borders of the surfaces with this surface class are \ntransparent to all surface molecules. \n   \n\n- ABSORPTIVE: If affected_complex_pattern refers to a volume molecule it is destroyed if it touches surfaces with this surface class. \nIf affected_complex_pattern refers to a surface molecule it is destroyed if it touches the border of the surface with \nthis surface class, i.e., it is allowed to release surface molecules on absorptive surfaces, they get destroyed only\nwhen they touch the border of this surface. \nTick marks on name allow destruction from only one side of the surface for volume molecules or selective destruction \nfor surface molecules on the surfaces’s border based on their orientation. \nTo make a surface with this surface class absorptive to all volume molecules, m.ALL_VOLUME_MOLECULES \ncan be used for affected_complex_pattern. \nTo make a border of the surface with this surface class absorptive to all surface molecules,\nm.ALL_SURFACE_MOLECULES can be used for name. \nUsing m.ALL_MOLECULES as affected_complex_pattern has the effect that surfaces with this surface\nclass are absorptive for all volume molecules and borders of the surfaces with this surface class \nare absorptive for all surface molecules.\n\n- CONCENTRATION_CLAMP: Clamps concentration at a surface by periodically releasing molecules that correspond\nto the wall being a transparent boundary to the area with given concentration, \nand by absorbing all molecules that hit this surface.  \nThe molecules matching affected_complex_pattern are destroyed if they touch the surface (as if they\nhad passed through), and new molecules are created at the surface, as if molecules had passed through \nfrom the other side at a concentration value (units = M). \nOrientation marks may be used; in this case, the other side of the surface is reflective. \nNote that this command is only used to set the effective concentration of a volume molecule at a surface; \nit is not valid to specify a surface molecule. \n\n- FLUX_CLAMP: Clamps flux at a surface by periodically releasing molecules that correspond\nto the wall being a transparent boundary to the area with given concentration. \nThe clamped surface reflects these molecules. \n\n")
    .value("UNSET", SurfacePropertyType::UNSET)
    .value("REACTIVE", SurfacePropertyType::REACTIVE)
    .value("REFLECTIVE", SurfacePropertyType::REFLECTIVE)
    .value("TRANSPARENT", SurfacePropertyType::TRANSPARENT)
    .value("ABSORPTIVE", SurfacePropertyType::ABSORPTIVE)
    .value("CONCENTRATION_CLAMP", SurfacePropertyType::CONCENTRATION_CLAMP)
    .value("FLUX_CLAMP", SurfacePropertyType::FLUX_CLAMP)
    .export_values();
  py::enum_<ExprNodeType>(m, "ExprNodeType", py::arithmetic(), "Used internally to represent expression trees.\n- UNSET\n\n- LEAF\n\n- ADD\n\n- SUB\n\n")
    .value("UNSET", ExprNodeType::UNSET)
    .value("LEAF", ExprNodeType::LEAF)
    .value("ADD", ExprNodeType::ADD)
    .value("SUB", ExprNodeType::SUB)
    .export_values();
  py::enum_<RegionNodeType>(m, "RegionNodeType", py::arithmetic(), "Used internally to represent region trees.\n- UNSET\n\n- LEAF_GEOMETRY_OBJECT\n\n- LEAF_SURFACE_REGION\n\n- UNION\n\n- DIFFERENCE\n\n- INTERSECT\n\n")
    .value("UNSET", RegionNodeType::UNSET)
    .value("LEAF_GEOMETRY_OBJECT", RegionNodeType::LEAF_GEOMETRY_OBJECT)
    .value("LEAF_SURFACE_REGION", RegionNodeType::LEAF_SURFACE_REGION)
    .value("UNION", RegionNodeType::UNION)
    .value("DIFFERENCE", RegionNodeType::DIFFERENCE)
    .value("INTERSECT", RegionNodeType::INTERSECT)
    .export_values();
  py::enum_<ReactionType>(m, "ReactionType", py::arithmetic(), "Used in reaction callbacks.\n- UNSET\n\n- UNIMOL_VOLUME\n\n- UNIMOL_SURFACE\n\n- VOLUME_VOLUME\n\n- VOLUME_SURFACE\n\n- SURFACE_SURFACE\n\n")
    .value("UNSET", ReactionType::UNSET)
    .value("UNIMOL_VOLUME", ReactionType::UNIMOL_VOLUME)
    .value("UNIMOL_SURFACE", ReactionType::UNIMOL_SURFACE)
    .value("VOLUME_VOLUME", ReactionType::VOLUME_VOLUME)
    .value("VOLUME_SURFACE", ReactionType::VOLUME_SURFACE)
    .value("SURFACE_SURFACE", ReactionType::SURFACE_SURFACE)
    .export_values();
  py::enum_<MoleculeType>(m, "MoleculeType", py::arithmetic(), "Used in molecule introspection and internally in checkpointing.\n- UNSET\n\n- VOLUME\n\n- SURFACE\n\n")
    .value("UNSET", MoleculeType::UNSET)
    .value("VOLUME", MoleculeType::VOLUME)
    .value("SURFACE", MoleculeType::SURFACE)
    .export_values();
}

} // namespace API
} // namespace MCell

