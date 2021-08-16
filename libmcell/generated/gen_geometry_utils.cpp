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

#include <sstream>
#include "api/pybind11_stl_include.h"
#include "api/python_export_utils.h"
#include "gen_geometry_utils.h"
#include "api/geometry_object.h"
#include "api/model.h"

namespace MCell {
namespace API {

void define_pybinding_geometry_utils(py::module& m) {
  m.def_submodule("geometry_utils")
      .def("create_box", &geometry_utils::create_box, py::arg("name"), py::arg("edge_dimension") = FLT_UNSET, py::arg("xyz_dimensions") = std::vector<double>(), "Creates a GeometryObject in the shape of a cube whose center is at (0, 0, 0).\n- name: Name of the created geometry object.\n\n- edge_dimension: Specifies length of each edge of the box in um. \nNone of x/y/z dimensions can be set.\n\n\n- xyz_dimensions: Specifies x/y/z sizes of the box in um. Parameter edge_dimension must not be set.\n\n")
      .def("create_icosphere", &geometry_utils::create_icosphere, py::arg("name"), py::arg("radius"), py::arg("subdivisions"), "Creates a GeometryObject in the shape of an icosphere whose center is at (0, 0, 0).\n- name: Name of the created geometry object.\n\n- radius: Specifies radius of the sphere.\n\n- subdivisions: Number of subdivisions from the initial icosphere. \nThe higher this value will be the smoother the icosphere will be.\nAllowed range is between 1 and 8.\n\n\n")
      .def("validate_volumetric_mesh", &geometry_utils::validate_volumetric_mesh, py::arg("model"), py::arg("geometry_object"), "Checks that the mesh was correctly analyzed, that it has volume and \nall edges have neighboring walls.\nMust be called after model initialization. \nThrows exception with detained message if validation did not pass. \n\n- model: Model object after initialization.\n\n- geometry_object: Geometry object to be checked.\n\n")
    ;
}

} // namespace API
} // namespace MCell

