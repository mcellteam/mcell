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
const std::string ALL_MOLECULES = "ALL_MOLECULES";
const std::string ALL_VOLUME_MOLECULES = "ALL_VOLUME_MOLECULES";
const std::string ALL_SURFACE_MOLECULES = "ALL_SURFACE_MOLECULES";

enum class Orientation {
  DOWN = -1,
  NONE = 0,
  UP = 1,
  NOT_SET = 2,
  ANY = 3
};


static inline  std::ostream& operator << (std::ostream& out, const Orientation v) {
  switch (v) {
    case Orientation::DOWN: out << "Orientation.DOWN (-1)"; break;
    case Orientation::NONE: out << "Orientation.NONE (0)"; break;
    case Orientation::UP: out << "Orientation.UP (1)"; break;
    case Orientation::NOT_SET: out << "Orientation.NOT_SET (2)"; break;
    case Orientation::ANY: out << "Orientation.ANY (3)"; break;
  }
  return out;
};


enum class Notification {
  NONE = 0,
  BRIEF = 1,
  FULL = 2
};


static inline  std::ostream& operator << (std::ostream& out, const Notification v) {
  switch (v) {
    case Notification::NONE: out << "Notification.NONE (0)"; break;
    case Notification::BRIEF: out << "Notification.BRIEF (1)"; break;
    case Notification::FULL: out << "Notification.FULL (2)"; break;
  }
  return out;
};


enum class WarningLevel {
  IGNORE = 0,
  WARNING = 1,
  ERROR = 2
};


static inline  std::ostream& operator << (std::ostream& out, const WarningLevel v) {
  switch (v) {
    case WarningLevel::IGNORE: out << "WarningLevel.IGNORE (0)"; break;
    case WarningLevel::WARNING: out << "WarningLevel.WARNING (1)"; break;
    case WarningLevel::ERROR: out << "WarningLevel.ERROR (2)"; break;
  }
  return out;
};


enum class VizMode {
  ASCII = 0,
  CELLBLENDER = 1
};


static inline  std::ostream& operator << (std::ostream& out, const VizMode v) {
  switch (v) {
    case VizMode::ASCII: out << "VizMode.ASCII (0)"; break;
    case VizMode::CELLBLENDER: out << "VizMode.CELLBLENDER (1)"; break;
  }
  return out;
};


enum class Shape {
  UNSET = 0,
  SPHERICAL = 1,
  REGION_EXPR = 2
};


static inline  std::ostream& operator << (std::ostream& out, const Shape v) {
  switch (v) {
    case Shape::UNSET: out << "Shape.UNSET (0)"; break;
    case Shape::SPHERICAL: out << "Shape.SPHERICAL (1)"; break;
    case Shape::REGION_EXPR: out << "Shape.REGION_EXPR (2)"; break;
  }
  return out;
};


enum class SurfacePropertyType {
  UNSET = 0,
  REFLECTIVE = 1,
  TRANSPARENT = 2,
  ABSORPTIVE = 3
};


static inline  std::ostream& operator << (std::ostream& out, const SurfacePropertyType v) {
  switch (v) {
    case SurfacePropertyType::UNSET: out << "SurfacePropertyType.UNSET (0)"; break;
    case SurfacePropertyType::REFLECTIVE: out << "SurfacePropertyType.REFLECTIVE (1)"; break;
    case SurfacePropertyType::TRANSPARENT: out << "SurfacePropertyType.TRANSPARENT (2)"; break;
    case SurfacePropertyType::ABSORPTIVE: out << "SurfacePropertyType.ABSORPTIVE (3)"; break;
  }
  return out;
};


enum class ExprNodeType {
  UNSET = 0,
  LEAF = 1,
  ADD = 2,
  SUB = 3
};


static inline  std::ostream& operator << (std::ostream& out, const ExprNodeType v) {
  switch (v) {
    case ExprNodeType::UNSET: out << "ExprNodeType.UNSET (0)"; break;
    case ExprNodeType::LEAF: out << "ExprNodeType.LEAF (1)"; break;
    case ExprNodeType::ADD: out << "ExprNodeType.ADD (2)"; break;
    case ExprNodeType::SUB: out << "ExprNodeType.SUB (3)"; break;
  }
  return out;
};


enum class RegionNodeType {
  UNSET = 0,
  LEAF_GEOMETRY_OBJECT = 1,
  LEAF_SURFACE_REGION = 2,
  UNION = 3,
  DIFFERENCE = 4,
  INTERSECT = 5
};


static inline  std::ostream& operator << (std::ostream& out, const RegionNodeType v) {
  switch (v) {
    case RegionNodeType::UNSET: out << "RegionNodeType.UNSET (0)"; break;
    case RegionNodeType::LEAF_GEOMETRY_OBJECT: out << "RegionNodeType.LEAF_GEOMETRY_OBJECT (1)"; break;
    case RegionNodeType::LEAF_SURFACE_REGION: out << "RegionNodeType.LEAF_SURFACE_REGION (2)"; break;
    case RegionNodeType::UNION: out << "RegionNodeType.UNION (3)"; break;
    case RegionNodeType::DIFFERENCE: out << "RegionNodeType.DIFFERENCE (4)"; break;
    case RegionNodeType::INTERSECT: out << "RegionNodeType.INTERSECT (5)"; break;
  }
  return out;
};


void define_pybinding_constants(py::module& m);
void define_pybinding_enums(py::module& m);

} // namespace API
} // namespace MCell

#endif // API_GEN_CONSTANTS

