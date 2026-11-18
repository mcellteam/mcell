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

#include "data_model_geometry.h"

#include <stdexcept>
#include <vector>

#include "geometry_analysis.hpp"
#include "mesh_table.hpp"

#include "include/datamodel_defines.h"
#include "generator_utils.h"

using namespace Json;
using Json::Value;

namespace MCell {

void compute_volume_and_area(Json::Value& model_object, double& volume, double& area) {

  string name = model_object[KEY_NAME].asString();

  Value& vertex_list = get_node(name, model_object, KEY_VERTEX_LIST);
  Value& element_connections = get_node(name, model_object, KEY_ELEMENT_CONNECTIONS);

  // the mesh in microns, as the data model defines it, assuming that the data
  // model is correct
  std::vector<GeometryAnalysis::Point> points;
  points.reserve(vertex_list.size());
  for (Value::ArrayIndex i = 0; i < vertex_list.size(); i++) {
    Value& vertex = vertex_list[i];
    points.push_back({vertex[0].asDouble(), vertex[1].asDouble(), vertex[2].asDouble()});
  }
  std::vector<GeometryAnalysis::Triangle> triangles;
  triangles.reserve(element_connections.size());
  for (Value::ArrayIndex i = 0; i < element_connections.size(); i++) {
    Value& element = element_connections[i];
    triangles.push_back({element[0].asInt(), element[1].asInt(), element[2].asInt()});
  }
  // the checks of MCell's converters on the mesh as given: walls are joined
  // only through their vertex indices, nothing is merged
  const MeshUtils::MeshCheck check = MeshUtils::check_mesh(points, triangles, points);
  if (!check.indices.out_of_range.empty()) {
    // the first index out of range, with the text this function gave before
    for (const int i: triangles[check.indices.out_of_range[0]]) {
      if (i < 0 || static_cast<std::size_t>(i) >= points.size()) {
        throw ConversionError("Triangle point index " + std::to_string(i) + " is out of range.");
      }
    }
  }
  if (!check.table_built) {
    // a wall that uses a vertex twice, or too many vertices or walls to check
    throw ConversionError("Geometry object " + name + " has invalid vertex indices, could not compute volume.");
  }
  // non-manifold edges and vertices are reported before edges of only one wall;
  // the orientation is not checked here, as before
  if (!check.manifold()) {
    throw ConversionError("Non-manifold mesh is not allowed, error for " + name + ".");
  }
  if (!MeshUtils::find_boundary_edges(check.table).empty()) {
    throw ConversionError("Geometry object " + name + " is not watertight, could not compute volume.");
  }

  // computed as MCell computes a region's volume and its walls' areas
  const GeometryAnalysis::VolumeArea volume_and_area = GeometryAnalysis::volume_and_area(points, triangles);
  volume = volume_and_area.volume;
  area = volume_and_area.area;
}

} // namespace MCell
