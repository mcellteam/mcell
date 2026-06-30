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
#include <nanobind/nanobind.h>
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

NB_MAKE_OPAQUE(std::vector<std::shared_ptr<MCell::API::BaseChkptMol>>)
NB_MAKE_OPAQUE(std::vector<std::shared_ptr<MCell::API::Complex>>)
NB_MAKE_OPAQUE(std::vector<std::shared_ptr<MCell::API::Component>>)
NB_MAKE_OPAQUE(std::vector<std::shared_ptr<MCell::API::ComponentType>>)
NB_MAKE_OPAQUE(std::vector<std::shared_ptr<MCell::API::Count>>)
NB_MAKE_OPAQUE(std::vector<std::shared_ptr<MCell::API::ElementaryMolecule>>)
NB_MAKE_OPAQUE(std::vector<std::shared_ptr<MCell::API::ElementaryMoleculeType>>)
NB_MAKE_OPAQUE(std::vector<std::shared_ptr<MCell::API::GeometryObject>>)
NB_MAKE_OPAQUE(std::vector<std::shared_ptr<MCell::API::InitialSurfaceRelease>>)
NB_MAKE_OPAQUE(std::vector<std::vector<double>>)
NB_MAKE_OPAQUE(std::vector<std::vector<int>>)
NB_MAKE_OPAQUE(std::vector<std::shared_ptr<MCell::API::MoleculeReleaseInfo>>)
NB_MAKE_OPAQUE(std::vector<std::shared_ptr<MCell::API::ReactionRule>>)
NB_MAKE_OPAQUE(std::vector<std::shared_ptr<MCell::API::ReleaseSite>>)
NB_MAKE_OPAQUE(std::vector<std::shared_ptr<MCell::API::Species>>)
NB_MAKE_OPAQUE(std::vector<std::shared_ptr<MCell::API::SurfaceClass>>)
NB_MAKE_OPAQUE(std::vector<std::shared_ptr<MCell::API::SurfaceProperty>>)
NB_MAKE_OPAQUE(std::vector<std::shared_ptr<MCell::API::SurfaceRegion>>)
NB_MAKE_OPAQUE(std::vector<std::shared_ptr<MCell::API::VizOutput>>)
NB_MAKE_OPAQUE(std::vector<double>)
NB_MAKE_OPAQUE(std::vector<int>)
NB_MAKE_OPAQUE(std::vector<std::string>)
NB_MAKE_OPAQUE(std::vector<uint64_t>)

#endif // GEN_VECTORS_MAKE_OPAQUE_H
