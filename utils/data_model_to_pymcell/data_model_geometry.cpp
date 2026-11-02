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

#include "geometry_analysis.h"

#include "include/datamodel_defines.h"
#include "generator_utils.h"

using namespace Json;
using Json::Value;

namespace MCell {

void compute_volume_and_area(Json::Value& model_object, double& volume, double& area) {

  string name = model_object[KEY_NAME].asString();

  Value& vertex_list = get_node(name, model_object, KEY_VERTEX_LIST);
  Value& element_connections = get_node(name, model_object, KEY_ELEMENT_CONNECTIONS);

  // the mesh in microns, as the data model defines it, cleaned (vertices with
  // equal coordinates merged), assuming that the data model is correct
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
  const GeometryAnalysis::Mesh mesh = GeometryAnalysis::clean_mesh(points, triangles);

  if (!GeometryAnalysis::is_watertight(mesh)) {
    throw ConversionError("Geometry object " + name + " is not watertight, could not compute volume.");
  }

  // computed as vtkPolyDataNormals and vtkMassProperties computed them
  const GeometryAnalysis::VolumeArea volume_and_area = GeometryAnalysis::volume_and_area_vtk(mesh);
  volume = volume_and_area.volume;
  area = volume_and_area.area;
}

} // namespace MCell
