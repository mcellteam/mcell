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
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<double>>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<int>>);
PYBIND11_MAKE_OPAQUE(std::vector<std::shared_ptr<MCell::API::MoleculeReleaseInfo>>);
PYBIND11_MAKE_OPAQUE(std::vector<std::shared_ptr<MCell::API::ReactionRule>>);
PYBIND11_MAKE_OPAQUE(std::vector<std::shared_ptr<MCell::API::ReleaseSite>>);
PYBIND11_MAKE_OPAQUE(std::vector<std::shared_ptr<MCell::API::Species>>);
PYBIND11_MAKE_OPAQUE(std::vector<std::shared_ptr<MCell::API::SurfaceClass>>);
PYBIND11_MAKE_OPAQUE(std::vector<std::shared_ptr<MCell::API::SurfaceProperty>>);
PYBIND11_MAKE_OPAQUE(std::vector<std::shared_ptr<MCell::API::SurfaceRegion>>);
PYBIND11_MAKE_OPAQUE(std::vector<std::shared_ptr<MCell::API::VizOutput>>);
PYBIND11_MAKE_OPAQUE(std::vector<double>);
PYBIND11_MAKE_OPAQUE(std::vector<int>);
PYBIND11_MAKE_OPAQUE(std::vector<std::string>);
PYBIND11_MAKE_OPAQUE(std::vector<uint64_t>);

#endif // GEN_VECTORS_MAKE_OPAQUE_H
