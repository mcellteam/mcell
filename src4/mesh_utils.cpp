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
//
// The geometry report (GeometryReport), whose lines come from the checks of
// each mesh (mesh_table.hpp), and the checks of how a model uses its meshes.
// Each problem these find gives one warning line
//   Warning: [CODE] <object name>: <text>
// which scripts can parse with ^Warning: \[([A-Z_]+)\] (\S+): (.*)$. The name is
// the geometry object's name as MCell has it (Scene.box for MDL). After a
// warning the run continues. The cases:
//
//   MESH_INCONSISTENT_ORIENTATION  objects whose volume a release by
//                                  concentration uses: two walls run along
//                                  their shared edge in the same direction
//   MESH_INSIDE_OUT                counted objects and compartments: a shell
//                                  faces the wrong way (see
//                                  GeometryAnalysis::judge_parts); a release
//                                  by concentration into an object whose
//                                  signed volume is negative
//   REGION_NOT_CLOSED              a release of volume molecules into a region
//                                  with an object that has an edge of only one
//                                  wall
//
// The facts come from each object's checks (GeometryObject::get_mesh_check),
// whose edges come from the vertex indices; the signed volume is recomputed
// after vertex moves. These checks run only when the checks of each mesh found
// no errors, so a mesh here is manifold and consistently oriented. For MDL
// models, MCell3's initialization also stops counted objects and release
// regions that are open or inconsistently oriented by its own rules first.

#include "mesh_utils.hpp"

#include <algorithm>
#include <cassert>
#include <unordered_map>

#include "logging.h"

#include "bng/bng_config.h"
#include "world.h"
#include "partition.h"
#include "geometry.h"
#include "release_event.h"
#include "region_expr.h"

namespace MCell {
namespace MeshUtils {

GeometryReport::GeometryReport(const std::string& file_name_, const double length_unit_)
  : file_name(file_name_), length_unit(length_unit_) {
  // the header; makes the directory if needed
  BNG::initialize_report_file(file_name, "Geometry");
  out.open(file_name, std::ios::out | std::ios::app);
}


void GeometryReport::add(const std::string& object_name, const MeshCheck& check) {
  if (out.is_open()) {
    write_report_block(out, object_name, check, length_unit);
  }
  for (const ProblemLine& line: problem_lines(check)) {
    (line.error ? error_lines : warning_lines).push_back(format_problem_line(line, object_name, file_name));
  }
  if (check.has_errors()) {
    objects_with_errors.push_back(object_name);
  }
}


std::string GeometryReport::finish() {
  out.close();
  for (const std::string& line: error_lines) {
    mcell_error_raw("%s\n", line.c_str());
  }
  for (const std::string& line: warning_lines) {
    mcell_error_raw("%s\n", line.c_str());
  }
  if (objects_with_errors.empty()) {
    return "";
  }
  return stop_message(objects_with_errors, file_name);
}


void get_mesh(
    const World* world, const GeometryObject& obj, const bool convert_units,
    std::vector<GeometryAnalysis::Point>& points, std::vector<GeometryAnalysis::Triangle>& triangles) {

  const Partition& p = world->get_partition(PARTITION_ID_INITIAL);

  points.clear();
  triangles.clear();
  triangles.reserve(obj.wall_indices.size());
  // only an index: the points are numbered in the order the walls first use them
  std::unordered_map<vertex_index_t, int> mesh_index;
  mesh_index.reserve(obj.wall_indices.size() * 2);
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
        it = mesh_index.emplace(w.vertex_indices[i], static_cast<int>(points.size())).first;
        points.push_back({pt.x, pt.y, pt.z});
      }
      t[i] = it->second;
    }
    triangles.push_back(t);
  }
}


GeometryAnalysis::Mesh get_unmerged_mesh(
    const World* world, const GeometryObject& obj, const bool convert_units) {

  std::vector<GeometryAnalysis::Point> points;
  std::vector<GeometryAnalysis::Triangle> triangles;
  get_mesh(world, obj, convert_units, points, triangles);
  return GeometryAnalysis::unmerged_mesh(points, triangles);
}


GeometryAnalysis::Mesh get_clean_mesh(
    const World* world, const GeometryObject& obj, const bool convert_units) {

  // the cleaned mesh depends only on the order of the walls
  std::vector<GeometryAnalysis::Point> points;
  std::vector<GeometryAnalysis::Triangle> triangles;
  get_mesh(world, obj, convert_units, points, triangles);
  return GeometryAnalysis::clean_mesh(points, triangles);
}


namespace {

double to_um3(const World* world, const double volume) {
  const double l = world->config.length_unit;
  return volume * l * l * l;
}


// the orientation warning, printed once per object
void warn_inconsistent_orientation(World* world, const GeometryObject& obj) {
  if (world->checked_meshes.orientation.insert(obj.id).second) {
    mcell_warn(
        "[MESH_INCONSISTENT_ORIENTATION] %s: One or more pairs of adjacent walls traverse(s) their shared "
        "edge(s) in the same direction, so their normals point to different sides of the surface.",
        obj.name.c_str());
  }
}


void check_geometry_object(World* world, const GeometryObject& obj) {
  if (!obj.is_counted_volume_or_compartment()) {
    return;
  }

  // watertight: initialize_counted_volumes stopped the run otherwise. The
  // shells that face the wrong way, and their signed volumes, from the object's
  // checks at initialization
  std::size_t shells = 0;
  double signed_volume = 0.0;
  for (const ShellFacts& f: obj.get_mesh_check().shell_facts) {
    if (f.wrong_facing) {
      shells++;
      signed_volume += f.signed_volume;
    }
  }
  if (shells == 0) {
    return;
  }
  if (obj.counted_volume_index_outside == COUNTED_VOLUME_INDEX_INTERSECTS) {
    // a molecule that crosses one of its walls gets its counted volume from its
    // position, not from the side of the wall it hits
    // (CollisionUtils::update_counted_volume_id_when_crossing_wall)
    mcell_warn(
        "[MESH_INSIDE_OUT] %s: the walls' normals point into the object (signed volume %g um^3).",
        obj.name.c_str(), to_um3(world, signed_volume));
  }
  else {
    mcell_warn(
        "[MESH_INSIDE_OUT] %s: the walls' normals point into the object (signed volume %g um^3), so molecules "
        "that enter it are taken as leaving and molecules that leave as entering.",
        obj.name.c_str(), to_um3(world, signed_volume));
  }
}


// the objects of a region expression, each once, in the order they appear
void collect_objects(const RegionExprNode* node, std::vector<geometry_object_id_t>& ids) {
  if (node == nullptr) {
    return;
  }
  if (node->op == RegionExprOperator::LEAF_GEOMETRY_OBJECT) {
    if (std::find(ids.begin(), ids.end(), node->geometry_object_id) == ids.end()) {
      ids.push_back(node->geometry_object_id);
    }
    return;
  }
  collect_objects(node->left, ids);
  collect_objects(node->right, ids);
}

void check_release(World* world, const ReleaseEvent& release) {
  if (release.release_shape != ReleaseShape::REGION || release.region_expr.root == nullptr ||
      world->get_all_species().get(release.species_id).is_surf()) {
    return; // not a release of volume molecules into a region
  }

  const RegionExprNode* root = release.region_expr.root;
  std::vector<geometry_object_id_t> ids;
  collect_objects(root, ids);

  // what the release does with its objects (ReleaseEvent::release_inside_regions):
  // a release by concentration into one object computes its number of molecules
  // from the object's volume; every release places its molecules with an inside
  // test, or removes molecules with one when its number is negative
  const bool by_concentration = release.release_number_method == ReleaseNumberMethod::CONCENTRATION_NUM;
  const bool uses_volume = by_concentration && root->op == RegionExprOperator::LEAF_GEOMETRY_OBJECT;
  const bool removes =
      (release.release_number_method == ReleaseNumberMethod::CONST_NUM && release.release_number < 0) ||
      (by_concentration && release.concentration < 0);
  const char* site = release.release_site_name.c_str();

  const Partition& p = world->get_partition(PARTITION_ID_INITIAL);
  for (const geometry_object_id_t id: ids) {
    if (!world->checked_meshes.release_regions.insert({release.release_site_name, id}).second) {
      continue;
    }
    const GeometryObject& obj = p.get_geometry_object_by_id(id);
    const MeshCheck& check = obj.get_mesh_check();

    if (!check.watertight()) {
      if (uses_volume) {
        mcell_warn(
            "[REGION_NOT_CLOSED] %s: release site %s computes its number of molecules from this object's volume "
            "and places them inside it, but the object is not closed, so neither its volume nor its inside has "
            "a meaning.",
            obj.name.c_str(), site);
      }
      else {
        mcell_warn(
            "[REGION_NOT_CLOSED] %s: release site %s %s molecules inside this object, but the object is not "
            "closed, so \"inside\" has no meaning.",
            obj.name.c_str(), site, removes ? "removes" : "places");
      }
      continue;
    }

    if (!uses_volume) {
      continue;
    }
    if (!check.consistently_oriented()) {
      warn_inconsistent_orientation(world, obj);
    }
    else {
      // the object's signed volume; until a vertex move it is the volume a
      // release by concentration uses, which Region computes with the same
      // function from the same walls in the same order
      const double volume = obj.get_signed_volume(p);
      if (volume < 0) {
        mcell_warn(
            "[MESH_INSIDE_OUT] %s: the walls' normals point into the object (signed volume %g um^3), so release "
            "site %s computes a negative number of molecules from its volume.",
            obj.name.c_str(), to_um3(world, volume), site);
      }
    }
  }
}

} // namespace


void check_release_regions(World* world, const ReleaseEvent& release) {
  check_release(world, release);
}


void check_meshes(World* world) {
  const Partition& p = world->get_partition(PARTITION_ID_INITIAL);
  for (const GeometryObject& obj: p.get_geometry_objects()) {
    check_geometry_object(world, obj);
  }

  std::vector<BaseEvent*> releases;
  world->scheduler.get_all_events_with_type_index(EVENT_TYPE_INDEX_RELEASE, releases);
  for (const BaseEvent* event: releases) {
    const ReleaseEvent* release = dynamic_cast<const ReleaseEvent*>(event);
    assert(release != nullptr);
    check_release(world, *release);
  }
}

} // namespace MeshUtils
} // namespace MCell
