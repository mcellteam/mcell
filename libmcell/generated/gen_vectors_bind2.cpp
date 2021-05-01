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

  py::bind_vector<std::vector<std::shared_ptr<MCell::API::Complex>>>(m,"VectorComplex");
  py::implicitly_convertible<py::list, std::vector<std::shared_ptr<MCell::API::Complex>>>();

  py::bind_vector<std::vector<std::shared_ptr<MCell::API::Component>>>(m,"VectorComponent");
  py::implicitly_convertible<py::list, std::vector<std::shared_ptr<MCell::API::Component>>>();

  py::bind_vector<std::vector<std::shared_ptr<MCell::API::ComponentType>>>(m,"VectorComponentType");
  py::implicitly_convertible<py::list, std::vector<std::shared_ptr<MCell::API::ComponentType>>>();

  py::bind_vector<std::vector<std::shared_ptr<MCell::API::Count>>>(m,"VectorCount");
  py::implicitly_convertible<py::list, std::vector<std::shared_ptr<MCell::API::Count>>>();

  py::bind_vector<std::vector<std::shared_ptr<MCell::API::ElementaryMolecule>>>(m,"VectorElementaryMolecule");
  py::implicitly_convertible<py::list, std::vector<std::shared_ptr<MCell::API::ElementaryMolecule>>>();

  py::bind_vector<std::vector<std::shared_ptr<MCell::API::ElementaryMoleculeType>>>(m,"VectorElementaryMoleculeType");
  py::implicitly_convertible<py::list, std::vector<std::shared_ptr<MCell::API::ElementaryMoleculeType>>>();

  py::bind_vector<std::vector<std::shared_ptr<MCell::API::GeometryObject>>>(m,"VectorGeometryObject");
  py::implicitly_convertible<py::list, std::vector<std::shared_ptr<MCell::API::GeometryObject>>>();

  py::bind_vector<std::vector<std::shared_ptr<MCell::API::InitialSurfaceRelease>>>(m,"VectorInitialSurfaceRelease");
  py::implicitly_convertible<py::list, std::vector<std::shared_ptr<MCell::API::InitialSurfaceRelease>>>();

  py::bind_vector<std::vector<std::vector<double>>>(m,"VectorVectorDouble");
  py::implicitly_convertible<py::list, std::vector<std::vector<double>>>();

  py::bind_vector<std::vector<std::vector<float>>>(m,"VectorVectorFloat");
  py::implicitly_convertible<py::list, std::vector<std::vector<float>>>();

  py::bind_vector<std::vector<std::vector<int>>>(m,"VectorVectorInt");
  py::implicitly_convertible<py::list, std::vector<std::vector<int>>>();

  py::bind_vector<std::vector<std::shared_ptr<MCell::API::MoleculeReleaseInfo>>>(m,"VectorMoleculeReleaseInfo");
  py::implicitly_convertible<py::list, std::vector<std::shared_ptr<MCell::API::MoleculeReleaseInfo>>>();

  py::bind_vector<std::vector<std::shared_ptr<MCell::API::ReactionRule>>>(m,"VectorReactionRule");
  py::implicitly_convertible<py::list, std::vector<std::shared_ptr<MCell::API::ReactionRule>>>();

  py::bind_vector<std::vector<std::shared_ptr<MCell::API::ReleaseSite>>>(m,"VectorReleaseSite");
  py::implicitly_convertible<py::list, std::vector<std::shared_ptr<MCell::API::ReleaseSite>>>();

  py::bind_vector<std::vector<std::shared_ptr<MCell::API::Species>>>(m,"VectorSpecies");
  py::implicitly_convertible<py::list, std::vector<std::shared_ptr<MCell::API::Species>>>();

  py::bind_vector<std::vector<std::shared_ptr<MCell::API::SurfaceClass>>>(m,"VectorSurfaceClass");
  py::implicitly_convertible<py::list, std::vector<std::shared_ptr<MCell::API::SurfaceClass>>>();

  py::bind_vector<std::vector<std::shared_ptr<MCell::API::SurfaceProperty>>>(m,"VectorSurfaceProperty");
  py::implicitly_convertible<py::list, std::vector<std::shared_ptr<MCell::API::SurfaceProperty>>>();

  py::bind_vector<std::vector<std::shared_ptr<MCell::API::SurfaceRegion>>>(m,"VectorSurfaceRegion");
  py::implicitly_convertible<py::list, std::vector<std::shared_ptr<MCell::API::SurfaceRegion>>>();

  py::bind_vector<std::vector<MCell::Vec3>>(m,"VectorVec3");
  py::implicitly_convertible<py::list, std::vector<MCell::Vec3>>();

  py::bind_vector<std::vector<std::shared_ptr<MCell::API::VizOutput>>>(m,"VectorVizOutput");
  py::implicitly_convertible<py::list, std::vector<std::shared_ptr<MCell::API::VizOutput>>>();

  py::bind_vector<std::vector<double>>(m,"VectorDouble");
  //py::implicitly_convertible<py::list, std::vector<double>>();

  py::bind_vector<std::vector<float>>(m,"VectorFloat");
  py::implicitly_convertible<py::list, std::vector<float>>();

  py::bind_vector<std::vector<int>>(m,"VectorInt");
  py::implicitly_convertible<py::list, std::vector<int>>();

  py::bind_vector<std::vector<std::string>>(m,"VectorStr");
  py::implicitly_convertible<py::list, std::vector<std::string>>();

  py::bind_vector<std::vector<uint64_t>>(m,"VectorUint64");
  py::implicitly_convertible<py::list, std::vector<uint64_t>>();

}

} // namespace API
} // namespace MCell

