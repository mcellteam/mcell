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
#include <climits>
#include <cfloat>
#include "api/globals.h"

namespace MCell {
namespace API {

const std::string STATE_UNSET = "STATE_UNSET";
const int STATE_UNSET_INT = -1;
const int BOND_UNBOUND = -1;
const int BOND_BOUND = -2;
const int BOND_ANY = -3;
const double PARTITION_EDGE_EXTRA_MARGIN_UM = 0.01;
const int DEFAULT_COUNT_BUFFER_SIZE = 100;
const std::string ALL_MOLECULES = "ALL_MOLECULES";
const std::string ALL_VOLUME_MOLECULES = "ALL_VOLUME_MOLECULES";
const std::string ALL_SURFACE_MOLECULES = "ALL_SURFACE_MOLECULES";
const std::string DEFAULT_CHECKPOINTS_DIR = "checkpoints";
const std::string DEFAULT_SEED_DIR_PREFIX = "seed_";
const int DEFAULT_SEED_DIR_DIGITS = 5;
const std::string DEFAULT_ITERATION_DIR_PREFIX = "it_";
const int ID_INVALID = -1;
const int NUMBER_OF_TRAINS_UNLIMITED = -1;
const double TIME_INFINITY = 1e140;
const int INT_UNSET = INT32_MAX;
const double FLT_UNSET = FLT_MAX;
const int RNG_SIZE = 256;

enum class Orientation {
  DOWN = -1,
  NONE = 0,
  UP = 1,
  NOT_SET = 2,
  ANY = 3,
  DEFAULT = 4
};


static inline std::ostream& operator << (std::ostream& out, const Orientation v) {
  switch (v) {
    case Orientation::DOWN: out << "m.Orientation.DOWN"; break;
    case Orientation::NONE: out << "m.Orientation.NONE"; break;
    case Orientation::UP: out << "m.Orientation.UP"; break;
    case Orientation::NOT_SET: out << "m.Orientation.NOT_SET"; break;
    case Orientation::ANY: out << "m.Orientation.ANY"; break;
    case Orientation::DEFAULT: out << "m.Orientation.DEFAULT"; break;
  }
  return out;
};


enum class Notification {
  NONE = 0,
  BRIEF = 1,
  FULL = 2
};


static inline std::ostream& operator << (std::ostream& out, const Notification v) {
  switch (v) {
    case Notification::NONE: out << "m.Notification.NONE"; break;
    case Notification::BRIEF: out << "m.Notification.BRIEF"; break;
    case Notification::FULL: out << "m.Notification.FULL"; break;
  }
  return out;
};


enum class WarningLevel {
  IGNORE = 0,
  WARNING = 1,
  ERROR = 2
};


static inline std::ostream& operator << (std::ostream& out, const WarningLevel v) {
  switch (v) {
    case WarningLevel::IGNORE: out << "m.WarningLevel.IGNORE"; break;
    case WarningLevel::WARNING: out << "m.WarningLevel.WARNING"; break;
    case WarningLevel::ERROR: out << "m.WarningLevel.ERROR"; break;
  }
  return out;
};


enum class VizMode {
  ASCII = 0,
  CELLBLENDER_V1 = 1,
  CELLBLENDER = 2
};


static inline std::ostream& operator << (std::ostream& out, const VizMode v) {
  switch (v) {
    case VizMode::ASCII: out << "m.VizMode.ASCII"; break;
    case VizMode::CELLBLENDER_V1: out << "m.VizMode.CELLBLENDER_V1"; break;
    case VizMode::CELLBLENDER: out << "m.VizMode.CELLBLENDER"; break;
  }
  return out;
};


enum class Shape {
  UNSET = 0,
  SPHERICAL = 1,
  REGION_EXPR = 2,
  LIST = 3,
  COMPARTMENT = 4
};


static inline std::ostream& operator << (std::ostream& out, const Shape v) {
  switch (v) {
    case Shape::UNSET: out << "m.Shape.UNSET"; break;
    case Shape::SPHERICAL: out << "m.Shape.SPHERICAL"; break;
    case Shape::REGION_EXPR: out << "m.Shape.REGION_EXPR"; break;
    case Shape::LIST: out << "m.Shape.LIST"; break;
    case Shape::COMPARTMENT: out << "m.Shape.COMPARTMENT"; break;
  }
  return out;
};


enum class SurfacePropertyType {
  UNSET = 0,
  REACTIVE = 1,
  REFLECTIVE = 2,
  TRANSPARENT = 3,
  ABSORPTIVE = 4,
  CONCENTRATION_CLAMP = 5,
  FLUX_CLAMP = 6
};


static inline std::ostream& operator << (std::ostream& out, const SurfacePropertyType v) {
  switch (v) {
    case SurfacePropertyType::UNSET: out << "m.SurfacePropertyType.UNSET"; break;
    case SurfacePropertyType::REACTIVE: out << "m.SurfacePropertyType.REACTIVE"; break;
    case SurfacePropertyType::REFLECTIVE: out << "m.SurfacePropertyType.REFLECTIVE"; break;
    case SurfacePropertyType::TRANSPARENT: out << "m.SurfacePropertyType.TRANSPARENT"; break;
    case SurfacePropertyType::ABSORPTIVE: out << "m.SurfacePropertyType.ABSORPTIVE"; break;
    case SurfacePropertyType::CONCENTRATION_CLAMP: out << "m.SurfacePropertyType.CONCENTRATION_CLAMP"; break;
    case SurfacePropertyType::FLUX_CLAMP: out << "m.SurfacePropertyType.FLUX_CLAMP"; break;
  }
  return out;
};


enum class ExprNodeType {
  UNSET = 0,
  LEAF = 1,
  ADD = 2,
  SUB = 3
};


static inline std::ostream& operator << (std::ostream& out, const ExprNodeType v) {
  switch (v) {
    case ExprNodeType::UNSET: out << "m.ExprNodeType.UNSET"; break;
    case ExprNodeType::LEAF: out << "m.ExprNodeType.LEAF"; break;
    case ExprNodeType::ADD: out << "m.ExprNodeType.ADD"; break;
    case ExprNodeType::SUB: out << "m.ExprNodeType.SUB"; break;
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


static inline std::ostream& operator << (std::ostream& out, const RegionNodeType v) {
  switch (v) {
    case RegionNodeType::UNSET: out << "m.RegionNodeType.UNSET"; break;
    case RegionNodeType::LEAF_GEOMETRY_OBJECT: out << "m.RegionNodeType.LEAF_GEOMETRY_OBJECT"; break;
    case RegionNodeType::LEAF_SURFACE_REGION: out << "m.RegionNodeType.LEAF_SURFACE_REGION"; break;
    case RegionNodeType::UNION: out << "m.RegionNodeType.UNION"; break;
    case RegionNodeType::DIFFERENCE: out << "m.RegionNodeType.DIFFERENCE"; break;
    case RegionNodeType::INTERSECT: out << "m.RegionNodeType.INTERSECT"; break;
  }
  return out;
};


enum class ReactionType {
  UNSET = 0,
  UNIMOL_VOLUME = 1,
  UNIMOL_SURFACE = 2,
  VOLUME_VOLUME = 3,
  VOLUME_SURFACE = 4,
  SURFACE_SURFACE = 5
};


static inline std::ostream& operator << (std::ostream& out, const ReactionType v) {
  switch (v) {
    case ReactionType::UNSET: out << "m.ReactionType.UNSET"; break;
    case ReactionType::UNIMOL_VOLUME: out << "m.ReactionType.UNIMOL_VOLUME"; break;
    case ReactionType::UNIMOL_SURFACE: out << "m.ReactionType.UNIMOL_SURFACE"; break;
    case ReactionType::VOLUME_VOLUME: out << "m.ReactionType.VOLUME_VOLUME"; break;
    case ReactionType::VOLUME_SURFACE: out << "m.ReactionType.VOLUME_SURFACE"; break;
    case ReactionType::SURFACE_SURFACE: out << "m.ReactionType.SURFACE_SURFACE"; break;
  }
  return out;
};


enum class MoleculeType {
  UNSET = 0,
  VOLUME = 1,
  SURFACE = 2
};


static inline std::ostream& operator << (std::ostream& out, const MoleculeType v) {
  switch (v) {
    case MoleculeType::UNSET: out << "m.MoleculeType.UNSET"; break;
    case MoleculeType::VOLUME: out << "m.MoleculeType.VOLUME"; break;
    case MoleculeType::SURFACE: out << "m.MoleculeType.SURFACE"; break;
  }
  return out;
};


enum class BNGSimulationMethod {
  NONE = 0,
  ODE = 1,
  SSA = 2,
  PLA = 3,
  NF = 4
};


static inline std::ostream& operator << (std::ostream& out, const BNGSimulationMethod v) {
  switch (v) {
    case BNGSimulationMethod::NONE: out << "m.BNGSimulationMethod.NONE"; break;
    case BNGSimulationMethod::ODE: out << "m.BNGSimulationMethod.ODE"; break;
    case BNGSimulationMethod::SSA: out << "m.BNGSimulationMethod.SSA"; break;
    case BNGSimulationMethod::PLA: out << "m.BNGSimulationMethod.PLA"; break;
    case BNGSimulationMethod::NF: out << "m.BNGSimulationMethod.NF"; break;
  }
  return out;
};


void define_pybinding_constants(py::module& m);
void define_pybinding_enums(py::module& m);

} // namespace API
} // namespace MCell

#endif // API_GEN_CONSTANTS

