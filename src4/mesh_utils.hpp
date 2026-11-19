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
// meshes (the usage errors), whose lines scripts can parse. The checks of each
// mesh are in mesh_table.hpp; the algorithms on plain meshes are in
// GeometryAnalysis; counted volumes are in CountedVolumeUtils.
//
// This header is valid C++17 so that the rest of MCell can include it; the
// implementation is C++23.

#include <fstream>
#include <string>
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

// the same mesh as given, nothing merged (GeometryAnalysis::unmerged_mesh), with
// its bounding box
GeometryAnalysis::Mesh get_unmerged_mesh(
    const World* world, const GeometryObject& obj, const bool convert_units = false);

// the same mesh, cleaned (vertices with equal coordinates merged)
GeometryAnalysis::Mesh get_clean_mesh(
    const World* world, const GeometryObject& obj, const bool convert_units = false);

// The usage errors of a model: uses of geometry objects that their meshes don't
// allow. Every such use concerns volume molecules, which have an inside only in
// a watertight object, and which find that inside from the way its walls face
// only when no shell faces the wrong way (MeshCheck::wrong_facing_shells). One
// line per object and use, e.g.
//   Error: [REGION_NOT_WATERTIGHT] Scene.box: the model counts inside this object, but it is non-watertight (4 boundary edges); see reports/geometry_report_00001.txt
//   Error: [MESH_INSIDE_OUT] Scene.box: the model counts inside this object, but 1 shell faces the wrong way; see reports/geometry_report_00001.txt
// finish prints every line, then the caller stops the run once with the message
// finish returns (usage_stop_message), which names every object with a usage
// error. The geometry report itself has no usage lines; the lines name it for
// the boundary edges and the shells.
class UsageErrors {
public:
  explicit UsageErrors(const std::string& report_file_name);

  // keeps the line of one use of the object if its mesh doesn't allow it.
  // use: what the model does with the object, e.g. USE_COUNTED, as the line of
  // a non-watertight object says it. facing_use: for a use that also needs the
  // object to face outward, what the line of a shell facing the wrong way says
  // the model does (for a release by concentration, which computes its number
  // from the object's signed volume, "<use> by concentration"); "" when the
  // facing doesn't matter
  void add(const std::string& object_name, const MeshCheck& check, const std::string& use,
      const std::string& facing_use);

  // prints the lines to MCell's error stream; returns the message of the stop,
  // or "" when there is no usage error
  std::string finish() const;

private:
  std::string report_file_name;
  std::vector<std::string> lines;
  std::vector<std::string> objects; // in the order of their first line
};

// the uses that aren't releases, as the lines name them
extern const char* const USE_COUNTED;     // "the model counts inside this object"
extern const char* const USE_COMPARTMENT; // "this object is a compartment"

// checks how the model uses its geometry objects: the counted objects (volume
// molecules and their reactions counted inside), the compartments, and the
// regions of every scheduled release of volume molecules. The objects come in
// the partition's order, which is the report's; for each, its counts, its
// compartment and then its releases in the scheduler's order. Prints the lines;
// returns the stop, or "" when every use is allowed. The converters call it once
// the releases and counts are converted; they stopped the run on mesh errors
// before, so every mesh here is manifold and consistently oriented.
std::string check_usage(World* world);

// the same for the region of one release, e.g. one that Model.release_molecules
// makes after the initialization
std::string check_release_usage(World* world, const ReleaseEvent& release);

} // namespace MeshUtils
} // namespace MCell

#endif // SRC4_MESH_UTILS_HPP_
