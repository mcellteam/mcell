/******************************************************************************
 *
 * Copyright (C) 2026 by
 * The Salk Institute for Biological Studies
 *
 * Use of this source code is governed by an MIT-style
 * license that can be found in the LICENSE file or at
 * https://opensource.org/licenses/MIT.
 *
******************************************************************************/

// Built as C++23 (see src4/CMakeLists.txt).

#include "mesh_utils.h"

#include <map>

#include "world.h"
#include "partition.h"
#include "geometry.h"

namespace MCell {
namespace MeshUtils {

void get_mesh(
    const World* world, const GeometryObject& obj, const bool convert_units,
    std::vector<GeometryAnalysis::Point>& points, std::vector<GeometryAnalysis::Triangle>& triangles) {

  const Partition& p = world->get_partition(PARTITION_ID_INITIAL);

  points.clear();
  triangles.clear();
  triangles.reserve(obj.wall_indices.size());
  std::map<vertex_index_t, int> mesh_index;
  for (wall_index_t wi: obj.wall_indices) {
    const Wall& w = p.get_wall(wi);
    GeometryAnalysis::Triangle t;
    for (uint i = 0; i < VERTICES_IN_TRIANGLE; i++) {
      auto it = mesh_index.find(w.vertex_indices[i]);
      if (it == mesh_index.end()) {
        Vec3 pt = p.get_geometry_vertex(w.vertex_indices[i]);
        if (convert_units) {
          pt = pt * Vec3(world->config.length_unit);
        }
        it = mesh_index.emplace(w.vertex_indices[i], (int)points.size()).first;
        points.push_back({pt.x, pt.y, pt.z});
      }
      t[i] = it->second;
    }
    triangles.push_back(t);
  }
}


GeometryAnalysis::Mesh get_clean_mesh(
    const World* world, const GeometryObject& obj, const bool convert_units) {

  // the cleaned mesh depends only on the order of the walls
  std::vector<GeometryAnalysis::Point> points;
  std::vector<GeometryAnalysis::Triangle> triangles;
  get_mesh(world, obj, convert_units, points, triangles);
  return GeometryAnalysis::clean_mesh(points, triangles);
}

} // namespace MeshUtils
} // namespace MCell
