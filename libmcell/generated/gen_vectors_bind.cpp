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


#include "api/pybind11_stl_include.h"
#include "generated/gen_vectors_make_opaque.h"

namespace py = nanobind;

#include "api/base_chkpt_mol.h"
#include "api/complex.h"
#include "api/component.h"
#include "api/component_type.h"
#include "api/count.h"
#include "api/elementary_molecule.h"
#include "api/elementary_molecule_type.h"
#include "api/geometry_object.h"
#include "api/initial_surface_release.h"
#include "api/molecule_release_info.h"
#include "api/reaction_rule.h"
#include "api/release_site.h"
#include "api/species.h"
#include "api/surface_class.h"
#include "api/surface_property.h"
#include "api/surface_region.h"
#include "api/viz_output.h"

namespace MCell {
namespace API {

void gen_vectors_bind(py::module_& m){
  py::bind_vector<std::vector<std::shared_ptr<MCell::API::BaseChkptMol>>>(m,"VectorBaseChkptMol");

  py::bind_vector<std::vector<std::shared_ptr<MCell::API::Complex>>>(m,"VectorComplex");

  py::bind_vector<std::vector<std::shared_ptr<MCell::API::Component>>>(m,"VectorComponent");

  py::bind_vector<std::vector<std::shared_ptr<MCell::API::ComponentType>>>(m,"VectorComponentType");

  py::bind_vector<std::vector<std::shared_ptr<MCell::API::Count>>>(m,"VectorCount");

  py::bind_vector<std::vector<std::shared_ptr<MCell::API::ElementaryMolecule>>>(m,"VectorElementaryMolecule");

  py::bind_vector<std::vector<std::shared_ptr<MCell::API::ElementaryMoleculeType>>>(m,"VectorElementaryMoleculeType");

  py::bind_vector<std::vector<std::shared_ptr<MCell::API::GeometryObject>>>(m,"VectorGeometryObject");

  py::bind_vector<std::vector<std::shared_ptr<MCell::API::InitialSurfaceRelease>>>(m,"VectorInitialSurfaceRelease");

  py::bind_vector<std::vector<std::vector<double>>>(m,"VectorVectorFloat");

  py::bind_vector<std::vector<std::vector<int>>>(m,"VectorVectorInt");

  py::bind_vector<std::vector<std::shared_ptr<MCell::API::MoleculeReleaseInfo>>>(m,"VectorMoleculeReleaseInfo");

  py::bind_vector<std::vector<std::shared_ptr<MCell::API::ReactionRule>>>(m,"VectorReactionRule");

  py::bind_vector<std::vector<std::shared_ptr<MCell::API::ReleaseSite>>>(m,"VectorReleaseSite");

  py::bind_vector<std::vector<std::shared_ptr<MCell::API::Species>>>(m,"VectorSpecies");

  py::bind_vector<std::vector<std::shared_ptr<MCell::API::SurfaceClass>>>(m,"VectorSurfaceClass");

  py::bind_vector<std::vector<std::shared_ptr<MCell::API::SurfaceProperty>>>(m,"VectorSurfaceProperty");

  py::bind_vector<std::vector<std::shared_ptr<MCell::API::SurfaceRegion>>>(m,"VectorSurfaceRegion");

  py::bind_vector<std::vector<std::shared_ptr<MCell::API::VizOutput>>>(m,"VectorVizOutput");

  py::bind_vector<std::vector<double>>(m,"VectorFloat");

  py::bind_vector<std::vector<int>>(m,"VectorInt");

  py::bind_vector<std::vector<std::string>>(m,"VectorStr");

  py::bind_vector<std::vector<uint64_t>>(m,"VectorUint64");

}

} // namespace API
} // namespace MCell

