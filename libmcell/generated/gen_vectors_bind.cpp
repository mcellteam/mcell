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
#include "pybind11/include/pybind11/stl_bind.h"

namespace py = pybind11;

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

void gen_vectors_bind(py::module& m){
  py::bind_vector<std::vector<std::shared_ptr<MCell::API::BaseChkptMol>>>(m,"VectorBaseChkptMol");
  py::implicitly_convertible<py::list, std::vector<std::shared_ptr<MCell::API::BaseChkptMol>>>();
  py::implicitly_convertible<py::tuple, std::vector<std::shared_ptr<MCell::API::BaseChkptMol>>>();

  py::bind_vector<std::vector<std::shared_ptr<MCell::API::Complex>>>(m,"VectorComplex");
  py::implicitly_convertible<py::list, std::vector<std::shared_ptr<MCell::API::Complex>>>();
  py::implicitly_convertible<py::tuple, std::vector<std::shared_ptr<MCell::API::Complex>>>();

  py::bind_vector<std::vector<std::shared_ptr<MCell::API::Component>>>(m,"VectorComponent");
  py::implicitly_convertible<py::list, std::vector<std::shared_ptr<MCell::API::Component>>>();
  py::implicitly_convertible<py::tuple, std::vector<std::shared_ptr<MCell::API::Component>>>();

  py::bind_vector<std::vector<std::shared_ptr<MCell::API::ComponentType>>>(m,"VectorComponentType");
  py::implicitly_convertible<py::list, std::vector<std::shared_ptr<MCell::API::ComponentType>>>();
  py::implicitly_convertible<py::tuple, std::vector<std::shared_ptr<MCell::API::ComponentType>>>();

  py::bind_vector<std::vector<std::shared_ptr<MCell::API::Count>>>(m,"VectorCount");
  py::implicitly_convertible<py::list, std::vector<std::shared_ptr<MCell::API::Count>>>();
  py::implicitly_convertible<py::tuple, std::vector<std::shared_ptr<MCell::API::Count>>>();

  py::bind_vector<std::vector<std::shared_ptr<MCell::API::ElementaryMolecule>>>(m,"VectorElementaryMolecule");
  py::implicitly_convertible<py::list, std::vector<std::shared_ptr<MCell::API::ElementaryMolecule>>>();
  py::implicitly_convertible<py::tuple, std::vector<std::shared_ptr<MCell::API::ElementaryMolecule>>>();

  py::bind_vector<std::vector<std::shared_ptr<MCell::API::ElementaryMoleculeType>>>(m,"VectorElementaryMoleculeType");
  py::implicitly_convertible<py::list, std::vector<std::shared_ptr<MCell::API::ElementaryMoleculeType>>>();
  py::implicitly_convertible<py::tuple, std::vector<std::shared_ptr<MCell::API::ElementaryMoleculeType>>>();

  py::bind_vector<std::vector<std::shared_ptr<MCell::API::GeometryObject>>>(m,"VectorGeometryObject");
  py::implicitly_convertible<py::list, std::vector<std::shared_ptr<MCell::API::GeometryObject>>>();
  py::implicitly_convertible<py::tuple, std::vector<std::shared_ptr<MCell::API::GeometryObject>>>();

  py::bind_vector<std::vector<std::shared_ptr<MCell::API::InitialSurfaceRelease>>>(m,"VectorInitialSurfaceRelease");
  py::implicitly_convertible<py::list, std::vector<std::shared_ptr<MCell::API::InitialSurfaceRelease>>>();
  py::implicitly_convertible<py::tuple, std::vector<std::shared_ptr<MCell::API::InitialSurfaceRelease>>>();

  py::bind_vector<std::vector<std::vector<double>>>(m,"VectorVectorFloat");
  py::implicitly_convertible<py::list, std::vector<std::vector<double>>>();
  py::implicitly_convertible<py::tuple, std::vector<std::vector<double>>>();

  py::bind_vector<std::vector<std::vector<int>>>(m,"VectorVectorInt");
  py::implicitly_convertible<py::list, std::vector<std::vector<int>>>();
  py::implicitly_convertible<py::tuple, std::vector<std::vector<int>>>();

  py::bind_vector<std::vector<std::shared_ptr<MCell::API::MoleculeReleaseInfo>>>(m,"VectorMoleculeReleaseInfo");
  py::implicitly_convertible<py::list, std::vector<std::shared_ptr<MCell::API::MoleculeReleaseInfo>>>();
  py::implicitly_convertible<py::tuple, std::vector<std::shared_ptr<MCell::API::MoleculeReleaseInfo>>>();

  py::bind_vector<std::vector<std::shared_ptr<MCell::API::ReactionRule>>>(m,"VectorReactionRule");
  py::implicitly_convertible<py::list, std::vector<std::shared_ptr<MCell::API::ReactionRule>>>();
  py::implicitly_convertible<py::tuple, std::vector<std::shared_ptr<MCell::API::ReactionRule>>>();

  py::bind_vector<std::vector<std::shared_ptr<MCell::API::ReleaseSite>>>(m,"VectorReleaseSite");
  py::implicitly_convertible<py::list, std::vector<std::shared_ptr<MCell::API::ReleaseSite>>>();
  py::implicitly_convertible<py::tuple, std::vector<std::shared_ptr<MCell::API::ReleaseSite>>>();

  py::bind_vector<std::vector<std::shared_ptr<MCell::API::Species>>>(m,"VectorSpecies");
  py::implicitly_convertible<py::list, std::vector<std::shared_ptr<MCell::API::Species>>>();
  py::implicitly_convertible<py::tuple, std::vector<std::shared_ptr<MCell::API::Species>>>();

  py::bind_vector<std::vector<std::shared_ptr<MCell::API::SurfaceClass>>>(m,"VectorSurfaceClass");
  py::implicitly_convertible<py::list, std::vector<std::shared_ptr<MCell::API::SurfaceClass>>>();
  py::implicitly_convertible<py::tuple, std::vector<std::shared_ptr<MCell::API::SurfaceClass>>>();

  py::bind_vector<std::vector<std::shared_ptr<MCell::API::SurfaceProperty>>>(m,"VectorSurfaceProperty");
  py::implicitly_convertible<py::list, std::vector<std::shared_ptr<MCell::API::SurfaceProperty>>>();
  py::implicitly_convertible<py::tuple, std::vector<std::shared_ptr<MCell::API::SurfaceProperty>>>();

  py::bind_vector<std::vector<std::shared_ptr<MCell::API::SurfaceRegion>>>(m,"VectorSurfaceRegion");
  py::implicitly_convertible<py::list, std::vector<std::shared_ptr<MCell::API::SurfaceRegion>>>();
  py::implicitly_convertible<py::tuple, std::vector<std::shared_ptr<MCell::API::SurfaceRegion>>>();

  py::bind_vector<std::vector<MCell::Vec3>>(m,"VectorVec3");
  py::implicitly_convertible<py::list, std::vector<MCell::Vec3>>();
  py::implicitly_convertible<py::tuple, std::vector<MCell::Vec3>>();

  py::bind_vector<std::vector<std::shared_ptr<MCell::API::VizOutput>>>(m,"VectorVizOutput");
  py::implicitly_convertible<py::list, std::vector<std::shared_ptr<MCell::API::VizOutput>>>();
  py::implicitly_convertible<py::tuple, std::vector<std::shared_ptr<MCell::API::VizOutput>>>();

  py::bind_vector<std::vector<double>>(m,"VectorFloat");
  py::implicitly_convertible<py::list, std::vector<double>>();
  py::implicitly_convertible<py::tuple, std::vector<double>>();

  py::bind_vector<std::vector<int>>(m,"VectorInt");
  py::implicitly_convertible<py::list, std::vector<int>>();
  py::implicitly_convertible<py::tuple, std::vector<int>>();

  py::bind_vector<std::vector<std::string>>(m,"VectorStr");
  py::implicitly_convertible<py::list, std::vector<std::string>>();
  py::implicitly_convertible<py::tuple, std::vector<std::string>>();

  py::bind_vector<std::vector<uint64_t>>(m,"VectorUint64");
  py::implicitly_convertible<py::list, std::vector<uint64_t>>();
  py::implicitly_convertible<py::tuple, std::vector<uint64_t>>();

}

} // namespace API
} // namespace MCell

