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

#ifndef API_GEN_CONSTANTS
#define API_GEN_CONSTANTS

#include <string>
#include "api/globals.h"

namespace MCell {
namespace API {

const std::string STATE_UNSET = "state_unset";
const int STATE_UNSET_INT = -1;
const int BOND_UNBOUND = 0;
const int BOND_BOUND = -1;
const float_t PARTITION_EDGE_EXTRA_MARGIN_UM = 0.01;
const int DEFAULT_COUNT_BUFFER_SIZE = 10000;
const std::string ALL_MOLECULES_NAME = "ALL_MOLECULES";
const std::string ALL_VOLUME_MOLECULES_NAME = "ALL_VOLUME_MOLECULES";
const std::string ALL_SURFACE_MOLECULES_NAME = "ALL_SURFACE_MOLECULES";

enum class Orientation {
  Down = -1,
  None = 0,
  Up = 1,
  NotSet = 2
};


static inline  std::ostream& operator << (std::ostream& out, const Orientation v) {
  switch (v) {
    case Orientation::Down: out << "Orientation.Down (-1)"; break;
    case Orientation::None: out << "Orientation.None (0)"; break;
    case Orientation::Up: out << "Orientation.Up (1)"; break;
    case Orientation::NotSet: out << "Orientation.NotSet (2)"; break;
  }
  return out;
};


enum class Notification {
  None = 0,
  Brief = 1,
  Full = 2
};


static inline  std::ostream& operator << (std::ostream& out, const Notification v) {
  switch (v) {
    case Notification::None: out << "Notification.None (0)"; break;
    case Notification::Brief: out << "Notification.Brief (1)"; break;
    case Notification::Full: out << "Notification.Full (2)"; break;
  }
  return out;
};


enum class WarningLevel {
  Ignore = 0,
  Warning = 1,
  Error = 2
};


static inline  std::ostream& operator << (std::ostream& out, const WarningLevel v) {
  switch (v) {
    case WarningLevel::Ignore: out << "WarningLevel.Ignore (0)"; break;
    case WarningLevel::Warning: out << "WarningLevel.Warning (1)"; break;
    case WarningLevel::Error: out << "WarningLevel.Error (2)"; break;
  }
  return out;
};


enum class VizMode {
  Ascii = 0,
  Cellblender = 1
};


static inline  std::ostream& operator << (std::ostream& out, const VizMode v) {
  switch (v) {
    case VizMode::Ascii: out << "VizMode.Ascii (0)"; break;
    case VizMode::Cellblender: out << "VizMode.Cellblender (1)"; break;
  }
  return out;
};


enum class Shape {
  Unset = 0,
  Spherical = 1,
  RegionExpr = 2
};


static inline  std::ostream& operator << (std::ostream& out, const Shape v) {
  switch (v) {
    case Shape::Unset: out << "Shape.Unset (0)"; break;
    case Shape::Spherical: out << "Shape.Spherical (1)"; break;
    case Shape::RegionExpr: out << "Shape.RegionExpr (2)"; break;
  }
  return out;
};


enum class ExprNodeType {
  Unset = 0,
  Leaf = 1,
  Add = 2,
  Sub = 3
};


static inline  std::ostream& operator << (std::ostream& out, const ExprNodeType v) {
  switch (v) {
    case ExprNodeType::Unset: out << "ExprNodeType.Unset (0)"; break;
    case ExprNodeType::Leaf: out << "ExprNodeType.Leaf (1)"; break;
    case ExprNodeType::Add: out << "ExprNodeType.Add (2)"; break;
    case ExprNodeType::Sub: out << "ExprNodeType.Sub (3)"; break;
  }
  return out;
};


enum class RegionNodeType {
  Unset = 0,
  LeafGeometryObject = 1,
  LeafSurfaceRegion = 2,
  Union = 3,
  Difference = 4,
  Intersect = 5
};


static inline  std::ostream& operator << (std::ostream& out, const RegionNodeType v) {
  switch (v) {
    case RegionNodeType::Unset: out << "RegionNodeType.Unset (0)"; break;
    case RegionNodeType::LeafGeometryObject: out << "RegionNodeType.LeafGeometryObject (1)"; break;
    case RegionNodeType::LeafSurfaceRegion: out << "RegionNodeType.LeafSurfaceRegion (2)"; break;
    case RegionNodeType::Union: out << "RegionNodeType.Union (3)"; break;
    case RegionNodeType::Difference: out << "RegionNodeType.Difference (4)"; break;
    case RegionNodeType::Intersect: out << "RegionNodeType.Intersect (5)"; break;
  }
  return out;
};


void define_pybinding_constants(py::module& m);
void define_pybinding_enums(py::module& m);

} // namespace API
} // namespace MCell

#endif // API_GEN_CONSTANTS

