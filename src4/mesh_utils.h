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

#ifndef SRC4_MESH_UTILS_H_
#define SRC4_MESH_UTILS_H_

// MCell's side of the mesh analysis: the meshes of geometry objects for
// GeometryAnalysis. The algorithms on plain meshes are in GeometryAnalysis;
// counted volumes are in CountedVolumeUtils.
//
// This header is valid C++17 so that the rest of MCell can include it; the
// implementation is C++23.

#include <vector>

#include "geometry_analysis.h"

namespace MCell {

class World;
class GeometryObject;

namespace MeshUtils {

// the walls of a geometry object in order, as triangles of the simulator's
// vertices in internal units (in microns with convert_units); each vertex is
// listed once, in the order the walls first use it
void get_mesh(
    const World* world, const GeometryObject& obj, const bool convert_units,
    std::vector<GeometryAnalysis::Point>& points, std::vector<GeometryAnalysis::Triangle>& triangles);

// the same mesh, cleaned (vertices with equal coordinates merged)
GeometryAnalysis::Mesh get_clean_mesh(
    const World* world, const GeometryObject& obj, const bool convert_units = false);

} // namespace MeshUtils
} // namespace MCell

#endif // SRC4_MESH_UTILS_H_
