/******************************************************************************
 *
 * Copyright (C) 2020 by
 * The Salk Institute for Biological Studies
 *
 * Use of this source code is governed by an MIT-style
 * license that can be found in the LICENSE file or at
 * https://opensource.org/licenses/MIT.
 *
******************************************************************************/

#include <memory>
#include <vector>

#include "compartment_utils.h"

#include "api/api_common.h"
#include "counted_volume_utils.h"
#include "geometry_object.h"
#include "geometry_analysis.hpp"
#include "mesh_table.hpp"

using namespace std;

namespace MCell {

using namespace CountedVolumeUtils;

namespace API {

// note: similar code is in MeshUtils::get_unmerged_mesh but uses different input data
static void convert_compartment_objects_to_geom_object_infos(
    const std::vector<std::shared_ptr<API::GeometryObject>>& compartment_objects,
    const double rcp_length_unit,
    const std::map<const API::GeometryObject*, MeshUtils::MeshCheck>& checks,
    GeomObjectInfoVector& compartment_infos) {

  for (const auto& obj_ptr: compartment_objects) {
    const API::GeometryObject& obj = *obj_ptr;

    // the object's flags, from its vertex indices (MeshUtils::check_lists): the
    // checks stopped the run on a mesh that isn't manifold or consistently
    // oriented, and MCell4Converter::check_compartments_watertight on a
    // compartment that isn't watertight
    const auto check_it = checks.find(&obj);
    release_assert(check_it != checks.end());
    release_assert(check_it->second.watertight());

    // the mesh as given, in internal units, with the coordinates the simulation
    // will use (the geometry objects are converted later); nothing is merged
    vector<GeometryAnalysis::Point> points;
    points.reserve(obj.vertex_list.size());
    for (size_t i = 0; i < obj.vertex_list.size(); i++) {
      assert(obj.vertex_list[i].size() == 3);
      const Vec3 v = obj.get_vertex_in_internal_units(i, rcp_length_unit);
      points.push_back({v.x, v.y, v.z});
    }
    vector<GeometryAnalysis::Triangle> triangles;
    triangles.reserve(obj.wall_list.size());
    for (const std::vector<int>& w: obj.wall_list) {
      assert(w.size() == VERTICES_IN_TRIANGLE);
      triangles.push_back({w[0], w[1], w[2]});
    }

    GeomObjectInfo obj_info(obj.name);
    obj_info.solid = make_shared<const GeometryAnalysis::Solid>(GeometryAnalysis::unmerged_mesh(points, triangles));
    compartment_infos.push_back(obj_info);
  }
}


static std::shared_ptr<API::GeometryObject> get_obj_by_name(
    std::vector<std::shared_ptr<API::GeometryObject>>& compartment_objects,
    const string& name) {
  for (auto& obj: compartment_objects) {
    if (obj->name == name) {
      return obj;
    }
  }
  release_assert(false);
  return std::shared_ptr<API::GeometryObject>(nullptr);
}


void set_parent_and_children_compartments(
    std::vector<std::shared_ptr<API::GeometryObject>>& compartment_objects,
    const double rcp_length_unit,
    const std::map<const API::GeometryObject*, MeshUtils::MeshCheck>& checks) {

  GeomObjectInfoVector compartment_infos;
  convert_compartment_objects_to_geom_object_infos(compartment_objects, rcp_length_unit, checks, compartment_infos);

  // get information on what objects are contained within each other
  // e.g. CP -> {EC} - CP is contained in EC, mapping contains all compartments, not just the direct parent
  ContainmentMap contained_in_mapping;
  IntersectingSet intersecting_object_infos;
  bool ok = compute_containement_mapping(
      nullptr, compartment_infos, contained_in_mapping, intersecting_object_infos);
  if (!ok) {
    throw RuntimeError("Unexpected error while computing compartment hierarchy.");
  }

  if (!intersecting_object_infos.empty()) {
    string names;
    for (const GeomObjectInfo& info: intersecting_object_infos) {
      names += info.name + ", ";
    }
    names = names.substr(0, names.size() - 2); // remove last comma
    throw RuntimeError("Error while computing compartment hierarchy, following compartments intersect " +
        names + ".");
  }

  // set children and parents
  PathLengthCache path_lengths;
  for (const auto& pair_name_parents: contained_in_mapping) {
    auto child = get_obj_by_name(compartment_objects, pair_name_parents.first.name);

    // select the direct parent
    const GeomObjectInfo* direct_parent_info =
        get_direct_parent_info(pair_name_parents.first, contained_in_mapping, path_lengths);

    if (direct_parent_info == nullptr) {
      // no parent, children set themselves as children of this compartment
    }
    else {
      auto parent = get_obj_by_name(compartment_objects, direct_parent_info->name);

      parent->child_compartments.insert(child);
      if (is_set(child->parent_compartment)) {
        // this should not really occur but lets better report it
        throw RuntimeError("Compartment object " + child->name + " cannot have multiple parent compartments " +
            child->parent_compartment->name + " and " + parent->name + ".");
      }
      child->parent_compartment = parent;
    }
  }
}


} // namespace API
} // namespace MCell
