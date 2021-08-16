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

#ifndef API_GEN_GEOMETRY_UTILS_H
#define API_GEN_GEOMETRY_UTILS_H

#include "api/api_common.h"

namespace MCell {
namespace API {

class GeometryObject;
class Model;
class PythonExportContext;

namespace geometry_utils {

std::shared_ptr<GeometryObject> create_box(const std::string& name, const double edge_dimension = FLT_UNSET, const std::vector<double> xyz_dimensions = std::vector<double>());
std::shared_ptr<GeometryObject> create_icosphere(const std::string& name, const double radius, const int subdivisions);
void validate_volumetric_mesh(std::shared_ptr<Model> model, std::shared_ptr<GeometryObject> geometry_object);

} // namespace geometry_utils

void define_pybinding_geometry_utils(py::module& m);
} // namespace API
} // namespace MCell

#endif // API_GEN_GEOMETRY_UTILS_H
