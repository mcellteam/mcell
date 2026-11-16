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

#ifndef SRC4_MESH_UTILS_HPP_
#define SRC4_MESH_UTILS_HPP_

// MCell's side of the mesh analysis: the geometry report, the meshes of
// geometry objects for GeometryAnalysis, and the checks of how a model uses its
// meshes, which print warnings that scripts can parse. The checks of each mesh
// are in mesh_table.hpp; the algorithms on plain meshes are in
// GeometryAnalysis; counted volumes are in CountedVolumeUtils.
//
// This header is valid C++17 so that the rest of MCell can include it; the
// implementation is C++23.

#include <fstream>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "defines.h"
#include "geometry_analysis.hpp"
#include "mesh_table.hpp"

namespace MCell {

class World;
class GeometryObject;
class ReleaseEvent;

namespace MeshUtils {

// The geometry report of a run, reports/geometry_report_<seed>.txt, and the
// terminal lines of the problems it lists. The converters check every geometry
// object before any wall is linked (check_lists for the Python API,
// check_parsed_walls for MDL) and add it here; the report gets one block per
// object (write_report_block) and is always written. finish then prints one line
// per object and problem to MCell's error stream, every error first, then every
// warning:
//   Error: [CODE] <object>: <text>; see reports/geometry_report_00001.txt
// The run goes on after warnings; after errors the caller stops it once, with
// the message finish returns, which names every object with errors.
class GeometryReport {
public:
  // starts the report file, with the header of BNG's reports
  GeometryReport(const std::string& file_name, const double length_unit);

  // writes the object's block and keeps its lines
  void add(const std::string& object_name, const MeshCheck& check);

  // prints the lines; returns the message of the stop, or "" when no object has
  // errors
  std::string finish();

private:
  std::string file_name;
  double length_unit;
  std::ofstream out;
  std::vector<std::string> error_lines;
  std::vector<std::string> warning_lines;
  std::vector<std::string> objects_with_errors;
};

// the walls of a geometry object in order, as triangles of the simulator's
// vertices in internal units (in microns with convert_units); each vertex is
// listed once, in the order the walls first use it
void get_mesh(
    const World* world, const GeometryObject& obj, const bool convert_units,
    std::vector<GeometryAnalysis::Point>& points, std::vector<GeometryAnalysis::Triangle>& triangles);

// the same mesh, cleaned (vertices with equal coordinates merged)
GeometryAnalysis::Mesh get_clean_mesh(
    const World* world, const GeometryObject& obj, const bool convert_units = false);

// what the checks have already looked at, so that each warning is printed
// once; World holds it
struct CheckedMeshes {
  // objects whose orientation was checked
  std::set<geometry_object_id_t> orientation;
  // release sites and the objects of their regions
  std::set<std::pair<std::string, geometry_object_id_t>> release_regions;
};

// checks how counted objects, compartments and the regions of all scheduled
// releases use their meshes, and prints one warning line for each problem found:
//   Warning: [CODE] <object name>: <text>
// The codes and the cases are listed in mesh_utils.cpp. The run continues.
// Must be called after the counted volumes are initialized.
void check_meshes(World* world);

// the checks of the regions of one release, e.g. of one made after the
// initialization
void check_release_regions(World* world, const ReleaseEvent& release);

} // namespace MeshUtils
} // namespace MCell

#endif // SRC4_MESH_UTILS_HPP_
