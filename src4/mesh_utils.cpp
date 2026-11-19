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
// each mesh (mesh_table.hpp), and the usage errors (UsageErrors): the uses of
// geometry objects that their meshes don't allow. Each usage error gives one
// line
//   Error: [CODE] <object name>: <use>, but <problem>; see <report>
// which scripts can parse with ^Error: \[([A-Z_]+)\] (\S+): (.*)$. The name is
// the geometry object's name as MCell has it (Scene.box for MDL). Every use
// concerns volume molecules; surface molecules are released onto walls and
// counted on walls, which needs no inside. The cases:
//
//   REGION_NOT_WATERTIGHT  the object is non-watertight (an edge of only one
//                          wall), and the model counts inside it, or it is a
//                          compartment, or a release of volume molecules
//                          places them inside it (or removes them from inside
//                          it) with the inside test (parity)
//   MESH_INSIDE_OUT        a shell of the object faces the wrong way (see
//                          GeometryAnalysis::judge_parts), and the model counts
//                          inside it (a molecule that crosses a wall takes the
//                          side the wall faces), or it is a compartment, or a
//                          release of volume molecules by concentration into
//                          the object alone computes its number from the
//                          object's signed volume
//
// The facts come from each object's checks (GeometryObject::get_mesh_check),
// whose edges come from the vertex indices. The converters check the uses once
// the releases and counts are converted, after they stopped the run on mesh
// errors, so a mesh here is manifold and consistently oriented. For MDL files,
// MCell3's initialization also stops counted objects and release regions that
// are open by its own rules first (it joins walls by their coordinates).

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


const char* const USE_COUNTED = "the model counts inside this object";
const char* const USE_COMPARTMENT = "this object is a compartment";


UsageErrors::UsageErrors(const std::string& report_file_name_)
  : report_file_name(report_file_name_) {
}


void UsageErrors::add(
    const std::string& object_name, const MeshCheck& check, const std::string& use, const std::string& facing_use) {

  // the checks of the meshes stopped the run on a mesh that isn't manifold or
  // consistently oriented, so the watertight check ran
  assert(check.watertight_checked);

  ProblemLine line;
  if (!check.watertight()) {
    const std::size_t n = check.boundary_edges;
    line.code = "REGION_NOT_WATERTIGHT";
    line.text = use + ", but it is non-watertight (" + std::to_string(n) +
        (n == 1 ? " boundary edge)" : " boundary edges)");
  }
  else if (!facing_use.empty() && check.wrong_facing_shells() != 0) {
    const std::size_t n = check.wrong_facing_shells();
    line.code = "MESH_INSIDE_OUT";
    line.text = facing_use + ", but " + std::to_string(n) + (n == 1 ? " shell faces" : " shells face") +
        " the wrong way";
  }
  else {
    return;
  }

  lines.push_back(format_problem_line(line, object_name, report_file_name));
  if (std::find(objects.begin(), objects.end(), object_name) == objects.end()) {
    objects.push_back(object_name);
  }
}


std::string UsageErrors::finish() const {
  for (const std::string& line: lines) {
    mcell_error_raw("%s\n", line.c_str());
  }
  if (objects.empty()) {
    return "";
  }
  return usage_stop_message(objects, report_file_name);
}


namespace {

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


// what a release of volume molecules does with the objects of its region, in
// the words of its lines
struct ReleaseUse {
  std::vector<geometry_object_id_t> objects; // each once, in the order the region names them
  std::string use;        // e.g. "release site rel releases volume molecules into this object"
  std::string facing_use; // "<use> by concentration", or "" when the facing doesn't matter
};


// false when the release isn't one of volume molecules into a region
bool get_release_use(const World* world, const ReleaseEvent& release, ReleaseUse& res) {
  if (release.release_shape != ReleaseShape::REGION || release.region_expr.root == nullptr ||
      world->get_all_species().get(release.species_id).is_surf()) {
    return false;
  }

  const RegionExprNode* root = release.region_expr.root;
  res.objects.clear();
  collect_objects(root, res.objects);

  // what the release does with its objects (ReleaseEvent::release_inside_regions):
  // it places its molecules with the inside test of each object, or removes
  // molecules with it when its number is negative; a release by concentration
  // into one object also computes its number from the object's signed volume,
  // so it needs the object to face outward too
  const bool by_concentration = release.release_number_method == ReleaseNumberMethod::CONCENTRATION_NUM;
  const bool removes =
      (release.release_number_method == ReleaseNumberMethod::CONST_NUM && release.release_number < 0) ||
      (by_concentration && release.concentration < 0);
  res.use = "release site " + release.release_site_name +
      (removes ? " removes volume molecules from this object" : " releases volume molecules into this object");
  const bool uses_volume = by_concentration && root->op == RegionExprOperator::LEAF_GEOMETRY_OBJECT;
  res.facing_use = uses_volume ? res.use + " by concentration" : "";
  return true;
}

} // namespace


std::string check_usage(World* world) {
  const Partition& p = world->get_partition(PARTITION_ID_INITIAL);

  // the releases of volume molecules into a region, in the scheduler's order
  std::vector<BaseEvent*> events;
  world->scheduler.get_all_events_with_type_index(EVENT_TYPE_INDEX_RELEASE, events);
  std::vector<ReleaseUse> releases;
  for (const BaseEvent* event: events) {
    const ReleaseEvent* release = dynamic_cast<const ReleaseEvent*>(event);
    assert(release != nullptr);
    ReleaseUse use;
    if (get_release_use(world, *release, use)) {
      releases.push_back(std::move(use));
    }
  }

  // each object's uses: its counts, its compartment, then its releases
  UsageErrors errors(world->config.get_geometry_report_file_name());
  for (const GeometryObject& obj: p.get_geometry_objects()) {
    const MeshCheck& check = obj.get_mesh_check();
    if (obj.is_counted_inside()) {
      errors.add(obj.name, check, USE_COUNTED, USE_COUNTED);
    }
    if (obj.represents_compartment()) {
      errors.add(obj.name, check, USE_COMPARTMENT, USE_COMPARTMENT);
    }
    for (const ReleaseUse& release: releases) {
      if (std::find(release.objects.begin(), release.objects.end(), obj.id) != release.objects.end()) {
        errors.add(obj.name, check, release.use, release.facing_use);
      }
    }
  }
  return errors.finish();
}


std::string check_release_usage(World* world, const ReleaseEvent& release) {
  UsageErrors errors(world->config.get_geometry_report_file_name());
  ReleaseUse use;
  if (get_release_use(world, release, use)) {
    const Partition& p = world->get_partition(PARTITION_ID_INITIAL);
    for (const geometry_object_id_t id: use.objects) {
      const GeometryObject& obj = p.get_geometry_object_by_id(id);
      errors.add(obj.name, obj.get_mesh_check(), use.use, use.facing_use);
    }
  }
  return errors.finish();
}

} // namespace MeshUtils
} // namespace MCell
