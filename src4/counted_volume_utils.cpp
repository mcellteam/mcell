/******************************************************************************
 *
 * Copyright (C) 2020 by
 * The Salk Institute for Biological Studies and
 * Pittsburgh Supercomputing Center, Carnegie Mellon University
 *
 * Use of this source code is governed by an MIT-style
 * license that can be found in the LICENSE file or at
 * https://opensource.org/licenses/MIT.
 *
******************************************************************************/

#include <algorithm>
#include <map>
#include <memory>

#include "logging.h"

#include "world.h"
#include "partition.h"
#include "geometry.h"

#include "counted_volume_utils.h"

using namespace std;

namespace MCell {

namespace CountedVolumeUtils {

GeometryObject& GeomObjectInfo::get_geometry_object_noconst(World* world) const {
  assert(world != nullptr);
  assert(for_counted_objects);
  Partition& p = world->get_partition(partition_id);
  return p.get_geometry_object(geometry_object_id);
}


const GeometryObject& GeomObjectInfo::get_geometry_object(const World* world) const {
  assert(world != nullptr);
  assert(for_counted_objects);
  const Partition& p = world->get_partition(partition_id);
  return p.get_geometry_object(geometry_object_id);
}


// the walls of a geometry object in order, as triangles of the simulator's
// vertices in internal units (in microns with convert_units); each vertex is
// listed once, in the order the walls first use it
static void get_mesh(
    const World* world, const GeometryObject& obj, const bool convert_units,
    vector<GeometryAnalysis::Point>& points, vector<GeometryAnalysis::Triangle>& triangles) {

  const Partition& p = world->get_partition(PARTITION_ID_INITIAL);

  points.clear();
  triangles.clear();
  triangles.reserve(obj.wall_indices.size());
  map<vertex_index_t, int> mesh_index;
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
  vector<GeometryAnalysis::Point> points;
  vector<GeometryAnalysis::Triangle> triangles;
  get_mesh(world, obj, convert_units, points, triangles);
  return GeometryAnalysis::clean_mesh(points, triangles);
}


// checks the mesh of every counted object and prepares it for the containment
// tests; reports every object that can't be counted in and then returns false
static bool prepare_counted_object_solids(const World* world, GeomObjectInfoVector& counted_objects) {
  bool res = true;

  for (GeomObjectInfo& obj_info: counted_objects) {

    const GeometryObject& obj = obj_info.get_geometry_object(world);
    GeometryAnalysis::Mesh mesh = get_clean_mesh(world, obj);

    // an edge shared by more than two walls is reported before an edge of only one wall
    if (GeometryAnalysis::has_non_manifold_edge(mesh)) {
      mcell_error_nodie("Non-manifold mesh is not allowed, error for %s.", obj.name.c_str());
      res = false;
    }
    else if (!GeometryAnalysis::is_watertight(mesh)) {
      mcell_error_nodie("Counting inside an open surface is not possible, error for %s.", obj.name.c_str());
      res = false;
    }
    else {
      obj_info.solid = make_shared<const GeometryAnalysis::Solid>(std::move(mesh));
    }
  }

  return res;
}


static const string& get_name(const World* world, const GeomObjectInfo& info) {
  if (world != nullptr) {
    return info.get_geometry_object(world).name;
  }
  else {
    return info.name;
  }
}


// returns false if theee was any error
bool compute_containement_mapping(
    const World* world, const GeomObjectInfoVector& counted_objects,
    ContainmentMap& contained_in_mapping,
    IntersectingSet& intersecting_objects
) {

  // keep it simple for now, let's just compute 'contained in' relation for each pair of objects
  // can be optimized in the future

  bool res = true;

  contained_in_mapping.clear();
  for (uint obj1 = 0; obj1 < counted_objects.size(); obj1++) {
    for (uint obj2 = obj1 + 1; obj2 < counted_objects.size(); obj2++) {
      assert(counted_objects[obj1].solid != nullptr && counted_objects[obj2].solid != nullptr);

      switch (GeometryAnalysis::containment(*counted_objects[obj1].solid, *counted_objects[obj2].solid)) {
        case GeometryAnalysis::Containment::FirstInSecond:
          contained_in_mapping[counted_objects[obj1]].insert(counted_objects[obj2]);
          break;

        case GeometryAnalysis::Containment::SecondInFirst:
          contained_in_mapping[counted_objects[obj2]].insert(counted_objects[obj1]);
          break;

        case GeometryAnalysis::Containment::Disjoint:
          // nothing to do
          break;

        case GeometryAnalysis::Containment::Intersect:
          // we are not processing all pairs so we need to insert both object ids
          intersecting_objects.insert(counted_objects[obj1]);
          intersecting_objects.insert(counted_objects[obj2]);
          break;

        case GeometryAnalysis::Containment::Identical:
          mcell_warn(
              "Identical counted objects are not supported yet, error for %s and %s.",
              get_name(world, counted_objects[obj1]).c_str(),
              get_name(world, counted_objects[obj2]).c_str()
          );
          res = false;
          break;
      }
    }
  }

  return res;
}


// length of the longest path from child through the objects that contain it to
// an object that no other object contains, the objects on the path counted
static int get_max_path_length(
    const GeomObjectInfo& child,
    const ContainmentMap& contained_in_mapping,
    PathLengthCache& path_lengths) {

  const auto cached = path_lengths.find(child);
  if (cached != path_lengths.end()) {
    return cached->second;
  }

  int max_length = 0;
  const auto it_obj_contained_in = contained_in_mapping.find(child);
  if (it_obj_contained_in != contained_in_mapping.end()) {
    for (const GeomObjectInfo& parent: it_obj_contained_in->second) {
      max_length = max(max_length, get_max_path_length(parent, contained_in_mapping, path_lengths));
    }
  }

  path_lengths[child] = max_length + 1;
  return max_length + 1;
}


const GeomObjectInfo* get_direct_parent_info(
    const GeomObjectInfo& obj_info, const ContainmentMap& contained_in_mapping,
    PathLengthCache& path_lengths) {

  // need to find an object that is a direct parent
  // we can find it by counting length of paths to from child to an object that has no parents
  // E.g., this is the contents of contained_in_mapping and we would like to determine parent of O3
  //    O3 -> O1, O2
  //    O2 -> O1
  //
  // Listing all possible paths from O3:
  //
  // 1: O3 -> O1
  // 2: O3 -> O2 -> O1
  //
  // the longest is 2:, therefore the direct parent of O3 is O2
  // partial overlaps and loops are not allowed here
  //
  // the lengths are cached, so each object's length is computed once;
  // with equal lengths, the first parent in the set's order wins

  const auto it_obj_contained_in = contained_in_mapping.find(obj_info);
  if (it_obj_contained_in == contained_in_mapping.end()) {
    return nullptr; // no parent
  }

  const GeomObjectInfo* best_parent = nullptr;
  int max_length = 0;
  for (const GeomObjectInfo& parent: it_obj_contained_in->second) {
    const int length = get_max_path_length(parent, contained_in_mapping, path_lengths);
    if (length > max_length) {
      max_length = length;
      best_parent = &parent;
    }
  }
  return best_parent;
}


static const GeometryObject* get_direct_parent(
    const World* world, const GeomObjectInfo& obj_info, const ContainmentMap& contained_in_mapping,
    PathLengthCache& path_lengths) {

  const GeomObjectInfo* info = get_direct_parent_info(obj_info, contained_in_mapping, path_lengths);
  if (info == nullptr) {
    return nullptr;
  }
  else {
    return &info->get_geometry_object(world);
  }
}


static counted_volume_index_t get_counted_volume_index_using_parents(
    World* world,
    const geometry_object_id_t obj_id,
    const ContainmentMap& contained_in_mapping,
    const IntersectingSet& intersecting_objects
) {
  // search in this map uses only object id
  GeomObjectInfo obj_info(PARTITION_ID_INVALID, obj_id);

  if (intersecting_objects.count(obj_info) != 0) {
    return COUNTED_VOLUME_INDEX_INTERSECTS;
  }

  const std::set<GeomObjectInfo>* obj_contained_in;

  auto it = contained_in_mapping.find(obj_info);
  if (it == contained_in_mapping.end()) {
    // not found - has no parent (intersection already checked)
    obj_contained_in = nullptr;
  }
  else {
    obj_contained_in = &it->second;
  }

  // we need to find or add a new counted volume to the partition
  Partition& p = world->get_partition(PARTITION_ID_INITIAL);

  // first prepare the counted volume
  CountedVolume cv;
  cv.contained_in_objects.insert(obj_id);
  if (obj_contained_in != nullptr) {
    for (const GeomObjectInfo& info: *obj_contained_in) {

      // if any of the parents intersects, let's assume that for counting this object intersects
      // as well, maybe some optimization can be done in the future
      if (intersecting_objects.count(info) != 0) {
        return COUNTED_VOLUME_INDEX_INTERSECTS;
      }

      geometry_object_index_t index = p.get_geometry_object(info.geometry_object_id).index;
      cv.contained_in_objects.insert_unique(index);
    }
  }

  return p.find_or_add_counted_volume(cv);
}


static void define_counted_volumes(
    World* world,
    const GeomObjectInfoVector& counted_objects,
    const ContainmentMap& contained_in_mapping,
    const IntersectingSet& intersecting_objects
) {
  PathLengthCache path_lengths;

  for (const GeomObjectInfo& obj_info: counted_objects) {

    GeometryObject& current_obj = obj_info.get_geometry_object_noconst(world);

    // set inside index for our object
    if (intersecting_objects.count(obj_info) != 0) {
      // intersects with other objects - not sure in which counted volume
      // a molecule ends up when entering this object
      current_obj.counted_volume_index_inside = COUNTED_VOLUME_INDEX_INTERSECTS;
    }
    else {
      // does not intersect - get counted volume from its parents
      current_obj.counted_volume_index_inside =
          get_counted_volume_index_using_parents(
              world, obj_info.geometry_object_id, contained_in_mapping, intersecting_objects);
    }

    const GeometryObject* direct_parent_obj =
        get_direct_parent(world, obj_info, contained_in_mapping, path_lengths);

    // set outside index for our object
    counted_volume_index_t outside_index;
    // parent was found
    if (direct_parent_obj != nullptr) {
      outside_index =
          get_counted_volume_index_using_parents(
              world, direct_parent_obj->id, contained_in_mapping, intersecting_objects);
    }
    else {
      // parent is unclear
      if (intersecting_objects.count(obj_info) == 0) {
        // this is a top level object that has no parent but still needs a counted volume to be created
        outside_index = COUNTED_VOLUME_INDEX_OUTSIDE_ALL;
      }
      else {
        // it is not clear what the outside is
        outside_index = COUNTED_VOLUME_INDEX_INTERSECTS;
      }
    }

    current_obj.counted_volume_index_outside = outside_index;
  }
}


// return true if counted volumes were correctly set up
// the only entry point to this utility file (so far)
bool initialize_counted_volumes(World* world, bool& has_intersecting_counted_objects) {

  GeomObjectInfoVector counted_objects;

  // get counted objects
  for (const Partition& p: world->get_partitions()) {
    for (const GeometryObject& obj: p.get_geometry_objects()) {
      if (obj.is_counted_volume_or_compartment()) {
        counted_objects.push_back( GeomObjectInfo(p.id, obj.id) );
      }
    }
  }

  // check the meshes and prepare them for the containment tests
  if (!prepare_counted_object_solids(world, counted_objects)) {
    return false;
  }

  // holds mapping that tells in which geometry objects is a given object contained
  ContainmentMap contained_in_mapping;

  // set of object ids that intersect
  IntersectingSet intersecting_objects;

  // compute 'contained-in' mapping
  bool ok = compute_containement_mapping(world, counted_objects, contained_in_mapping, intersecting_objects);
  if (!ok) {
    return false;
  }

  // define counted volumes, sets inside and outside volume id for geometry objects
  define_counted_volumes(world, counted_objects, contained_in_mapping, intersecting_objects);

  has_intersecting_counted_objects = !intersecting_objects.empty();

  return true;
}


// auxiliary function to compute volume, not related to counted volumes,
// returns FLT_INVALID if the object is not closed
double get_geometry_object_volume(const World* world, const GeometryObject& obj) {
  vector<GeometryAnalysis::Point> points;
  vector<GeometryAnalysis::Triangle> triangles;
  get_mesh(world, obj, false, points, triangles);
  if (!GeometryAnalysis::is_closed(GeometryAnalysis::clean_mesh(points, triangles))) {
    return FLT_INVALID;
  }

  // the volume of the object's region, which releases by concentration use:
  // the same function on the same walls in the same order
  return GeometryAnalysis::volume_and_area(points, triangles).volume;
}


void export_geometry_objects_to_obj(
    const World* world, const GeometryObjectVector& objs, const std::string& file_prefix) {

  // one object per geometry object, named by it, in microns; the material of
  // each is white with 25% transparency (note: same as DEFAULT_COLOR)
  vector<GeometryAnalysis::NamedMesh> meshes;
  for (const GeometryObject& o: objs) {
    meshes.push_back({o.name, get_clean_mesh(world, o, true)});
  }

  if (!GeometryAnalysis::write_obj(file_prefix, meshes)) {
    mcell_error_nodie("Could not write geometry to %s.obj and %s.mtl.", file_prefix.c_str(), file_prefix.c_str());
  }
}


} // namespace CountedVolumeUtils
} // namespace MCell
