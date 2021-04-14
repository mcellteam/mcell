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

#ifndef GEN_VECTORS_MAKE_OPAQUE_H
#define GEN_VECTORS_MAKE_OPAQUE_H

#include <vector>
#include <memory>
#include "pybind11/include/pybind11/pybind11.h"
#include "defines.h"

namespace MCell {
namespace API {

class BaseChkptMol;
class Complex;
class Component;
class ComponentType;
class Count;
class ElementaryMolecule;
class ElementaryMoleculeType;
class GeometryObject;
class InitialSurfaceRelease;
class MoleculeReleaseInfo;
class ReactionRule;
class ReleaseSite;
class Species;
class SurfaceClass;
class SurfaceProperty;
class SurfaceRegion;
class VizOutput;

} // namespace API
} // namespace MCell

PYBIND11_MAKE_OPAQUE(std::vector<std::shared_ptr<MCell::API::BaseChkptMol>>);
PYBIND11_MAKE_OPAQUE(std::vector<std::shared_ptr<MCell::API::Complex>>);
PYBIND11_MAKE_OPAQUE(std::vector<std::shared_ptr<MCell::API::Component>>);
PYBIND11_MAKE_OPAQUE(std::vector<std::shared_ptr<MCell::API::ComponentType>>);
PYBIND11_MAKE_OPAQUE(std::vector<std::shared_ptr<MCell::API::Count>>);
PYBIND11_MAKE_OPAQUE(std::vector<std::shared_ptr<MCell::API::ElementaryMolecule>>);
PYBIND11_MAKE_OPAQUE(std::vector<std::shared_ptr<MCell::API::ElementaryMoleculeType>>);
PYBIND11_MAKE_OPAQUE(std::vector<std::shared_ptr<MCell::API::GeometryObject>>);
PYBIND11_MAKE_OPAQUE(std::vector<std::shared_ptr<MCell::API::InitialSurfaceRelease>>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<float_t>>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<int>>);
PYBIND11_MAKE_OPAQUE(std::vector<std::shared_ptr<MCell::API::MoleculeReleaseInfo>>);
PYBIND11_MAKE_OPAQUE(std::vector<std::shared_ptr<MCell::API::ReactionRule>>);
PYBIND11_MAKE_OPAQUE(std::vector<std::shared_ptr<MCell::API::ReleaseSite>>);
PYBIND11_MAKE_OPAQUE(std::vector<std::shared_ptr<MCell::API::Species>>);
PYBIND11_MAKE_OPAQUE(std::vector<std::shared_ptr<MCell::API::SurfaceClass>>);
PYBIND11_MAKE_OPAQUE(std::vector<std::shared_ptr<MCell::API::SurfaceProperty>>);
PYBIND11_MAKE_OPAQUE(std::vector<std::shared_ptr<MCell::API::SurfaceRegion>>);
PYBIND11_MAKE_OPAQUE(std::vector<MCell::Vec3>);
PYBIND11_MAKE_OPAQUE(std::vector<std::shared_ptr<MCell::API::VizOutput>>);
PYBIND11_MAKE_OPAQUE(std::vector<float_t>);
PYBIND11_MAKE_OPAQUE(std::vector<int>);
PYBIND11_MAKE_OPAQUE(std::vector<std::string>);
PYBIND11_MAKE_OPAQUE(std::vector<uint64_t>);

#endif // GEN_VECTORS_MAKE_OPAQUE_H
