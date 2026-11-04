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
// The mesh checks. Each problem found gives one warning line
//   Warning: [CODE] <object name>: <text>
// which scripts can parse with ^Warning: \[([A-Z_]+)\] (\S+): (.*)$. The name is
// the geometry object's name as MCell has it (Scene.box for MDL). After a
// warning the run continues. The cases:
//
//   MESH_DUPLICATE_VERTICES        every object: two or more of the vertices
//                                  that its walls use have exactly the same
//                                  coordinates
//   MESH_INCONSISTENT_ORIENTATION  counted objects, compartments, and objects
//                                  whose volume a release by concentration
//                                  uses: two walls run along their shared edge
//                                  in the same direction
//   MESH_INSIDE_OUT                counted objects and compartments: a part
//                                  faces the wrong way (see
//                                  GeometryAnalysis::find_inside_out_parts); a
//                                  release by concentration into an object
//                                  whose signed volume is negative
//   REGION_NOT_CLOSED              a release of volume molecules into a region
//                                  with an object that has an edge of only one
//                                  wall
//
// Edges are matched by exact coordinates, as everywhere in GeometryAnalysis.
// Objects with an edge of three or more walls don't reach these checks in the
// Python API, which rejects them. MDL models are checked by MCell3 first: it
// stops on counted objects and release regions that are open or inconsistently
// oriented, so for MDL only MESH_DUPLICATE_VERTICES and MESH_INSIDE_OUT can
// appear.

#include "mesh_utils.h"

#include <algorithm>
#include <cassert>
#include <map>
#include <unordered_map>

#include "logging.h"

#include "world.h"
#include "partition.h"
#include "geometry.h"
#include "release_event.h"
#include "region_expr.h"

namespace MCell {
namespace MeshUtils {

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
  std::vector<GeometryAnalysis::Point> points;
  std::vector<GeometryAnalysis::Triangle> triangles;
  get_mesh(world, obj, false, points, triangles);

  if (GeometryAnalysis::count_duplicate_points(points) != 0) {
    mcell_warn(
        "[MESH_DUPLICATE_VERTICES] %s: Two or more vertices have exactly the same coordinates; MCell joins "
        "walls by coordinates, so it treats such vertices as one vertex.",
        obj.name.c_str());
  }

  if (!obj.is_counted_volume_or_compartment()) {
    return;
  }

  // closed and without an edge of three or more walls: initialize_counted_volumes
  // stopped the run otherwise
  const GeometryAnalysis::Mesh mesh = GeometryAnalysis::clean_mesh(points, triangles);
  if (GeometryAnalysis::count_inconsistent_edges(mesh) != 0) {
    warn_inconsistent_orientation(world, obj);
  }

  const GeometryAnalysis::InsideOutParts parts = GeometryAnalysis::find_inside_out_parts(mesh);
  if (parts.components == 0) {
    return;
  }
  if (obj.counted_volume_index_outside == COUNTED_VOLUME_INDEX_INTERSECTS) {
    // a molecule that crosses one of its walls gets its counted volume from its
    // position, not from the side of the wall it hits
    // (CollisionUtils::update_counted_volume_id_when_crossing_wall)
    mcell_warn(
        "[MESH_INSIDE_OUT] %s: the walls' normals point into the object (signed volume %g um^3).",
        obj.name.c_str(), to_um3(world, parts.signed_volume));
  }
  else {
    mcell_warn(
        "[MESH_INSIDE_OUT] %s: the walls' normals point into the object (signed volume %g um^3), so molecules "
        "that enter it are taken as leaving and molecules that leave as entering.",
        obj.name.c_str(), to_um3(world, parts.signed_volume));
  }
}


// what the release checks need to know about an object's mesh
struct RegionFacts {
  bool watertight;  // no edge of only one wall
  bool consistent;  // no two walls run along their shared edge in the same direction
  double volume;    // signed, as Region::initialize_volume_info_if_needed computes it
};

// one object's facts, computed once per cache; releases often share objects
// (in mdl/3640, 20 release sites use the same four objects of 369k walls)
const RegionFacts& get_region_facts(
    const World* world, const GeometryObject& obj, std::map<geometry_object_id_t, RegionFacts>& cache) {

  const auto it = cache.find(obj.id);
  if (it != cache.end()) {
    return it->second;
  }
  std::vector<GeometryAnalysis::Point> points;
  std::vector<GeometryAnalysis::Triangle> triangles;
  get_mesh(world, obj, false, points, triangles);
  const GeometryAnalysis::Mesh mesh = GeometryAnalysis::clean_mesh(points, triangles);
  // the volume that a release by concentration uses: Region computes it with the
  // same function from the same walls in the same order
  const RegionFacts facts{
      GeometryAnalysis::is_watertight(mesh),
      GeometryAnalysis::count_inconsistent_edges(mesh) == 0,
      GeometryAnalysis::volume_and_area(points, triangles).volume};
  return cache.emplace(obj.id, facts).first->second;
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

void check_release(
    World* world, const ReleaseEvent& release, std::map<geometry_object_id_t, RegionFacts>& cache) {
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
      (release.release_number_method == ReleaseNumberMethod::CONST_NUM && (int)release.release_number < 0) ||
      (by_concentration && release.concentration < 0);
  const char* site = release.release_site_name.c_str();

  const Partition& p = world->get_partition(PARTITION_ID_INITIAL);
  for (const geometry_object_id_t id: ids) {
    if (!world->checked_meshes.release_regions.insert({release.release_site_name, id}).second) {
      continue;
    }
    const GeometryObject& obj = p.get_geometry_object_by_id(id);
    const RegionFacts& facts = get_region_facts(world, obj, cache);

    if (!facts.watertight) {
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
    if (!facts.consistent) {
      warn_inconsistent_orientation(world, obj);
    }
    else if (facts.volume < 0) {
      mcell_warn(
          "[MESH_INSIDE_OUT] %s: the walls' normals point into the object (signed volume %g um^3), so release "
          "site %s computes a negative number of molecules from its volume.",
          obj.name.c_str(), to_um3(world, facts.volume), site);
    }
  }
}

} // namespace


void check_release_regions(World* world, const ReleaseEvent& release) {
  std::map<geometry_object_id_t, RegionFacts> cache;
  check_release(world, release, cache);
}


void check_meshes(World* world) {
  const Partition& p = world->get_partition(PARTITION_ID_INITIAL);
  for (const GeometryObject& obj: p.get_geometry_objects()) {
    check_geometry_object(world, obj);
  }

  std::vector<BaseEvent*> releases;
  world->scheduler.get_all_events_with_type_index(EVENT_TYPE_INDEX_RELEASE, releases);
  std::map<geometry_object_id_t, RegionFacts> cache;
  for (const BaseEvent* event: releases) {
    const ReleaseEvent* release = dynamic_cast<const ReleaseEvent*>(event);
    assert(release != nullptr);
    check_release(world, *release, cache);
  }
}

} // namespace MeshUtils
} // namespace MCell
