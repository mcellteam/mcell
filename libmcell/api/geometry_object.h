/******************************************************************************
 *
 * Copyright (C) 2020 by
 * The Salk Institute for Biological Studies
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
 * USA.
 *
******************************************************************************/

#ifndef API_GEOMETRY_OBJECT_H
#define API_GEOMETRY_OBJECT_H

#include "bng/bng_defines.h"

#include "generated/gen_geometry_object.h"
#include "api/common.h"
#include "api/surface_region.h"
#include "defines.h"

namespace MCell {
namespace API {

class GeometryObject;
typedef std::set<std::shared_ptr<API::GeometryObject>> GeometryObjectSet;

class GeometryObject: public GenGeometryObject {
public:
  GEOMETRY_OBJECT_CTOR()

public:
  void postprocess_in_ctor() override {
    node_type = RegionNodeType::LEAF_GEOMETRY_OBJECT;
    partition_id = PARTITION_ID_INVALID;
    geometry_object_id = GEOMETRY_OBJECT_ID_INVALID;
    first_vertex_index = VERTEX_INDEX_INVALID;
    parent_compartment = nullptr;
    vol_compartment_id = BNG::COMPARTMENT_ID_INVALID;
    surf_compartment_id = BNG::COMPARTMENT_ID_INVALID;

    for (auto& sr: surface_regions) {
      // not using shared pointers here, any attempt so far resulted in bad_weak_ptr exception
      // this is safe because the geometry object (parent) has a reference to the surface region
      sr->parent = this;
    }
  }

  // TODO: move to c++
  void check_semantics() const override {
    GenGeometryObject::check_semantics();
    for (auto& v: vertex_list) {
      if (v.size() != 3) {
        throw ValueError(
            S("Each item in the '") + NAME_VERTEX_LIST + "' argument must be a triplet of floats, error for " +
            vec_nonptr_to_str(v) + ".");
      }
    }

    for (auto& e: wall_list) {
      if (e.size() != 3) {
        throw ValueError(
            S("Each item in the '") + NAME_WALL_LIST + "' argument must be a triplet of integers, error for " +
            vec_nonptr_to_str(e) + ".");
        for (int vertex_index: e) {
          if (vertex_index < 0 || vertex_index >= (int)vertex_list.size()) {
            throw ValueError(
                S("Vertex index the '") + NAME_WALL_LIST + "' is out of range, error for " +
                std::to_string(vertex_index));
          }
        }
      }
    }

    for (auto& sr: surface_regions) {
      for (int wall_index: sr->wall_indices) {
        if (wall_index >= (int)wall_list.size()) {
          throw ValueError(
              S("Wall index in the '") + NAME_WALL_INDICES + "' of '" + sr->name + "' is out of range, error for " +
              std::to_string(wall_index));
        }
      }
    }
  }

  void translate(const Vec3& move) override {
    if (geometry_object_id != GEOMETRY_OBJECT_ID_INVALID) {
      throw RuntimeError(S("Method ") + NAME_TRANSLATE + " may be called only before model initialization.");
    }
    for (auto& v: vertex_list) {
      v[0] += move.x;
      v[1] += move.y;
      v[2] += move.z;
    }
  }

  // --- added manually ---
  void check_is_initialized() {
    if (geometry_object_id == GEOMETRY_OBJECT_ID_INVALID) {
      throw RuntimeError("Geometry object " + name + " is not present in model (or model was not initialized).");
    }
  }

  vertex_index_t get_partition_vertex_index(const int vertex_index) const {
    check_vertex_index(vertex_index);
    return first_vertex_index + vertex_index;
  }

  wall_index_t get_partition_wall_index(const int wall_index) const {
    check_wall_index(wall_index);
    return first_wall_index + wall_index;
  }

  int get_object_wall_index(const wall_index_t wall_index) const {
    int res = wall_index - first_wall_index;
    assert(res >= 0 && res < (int)wall_list.size());
    return res;
  }

private:
  void check_vertex_index(const int vertex_index) const {
    if (vertex_index < 0 || vertex_index >= (int)vertex_list.size()) {
      throw RuntimeError(
          "Vertex index " + std::to_string(vertex_index) + " is out of range for " + NAME_VERTEX_LIST + " of " + name + ".");
    }
  }

  void check_wall_index(const int wall_index) const {
    if (wall_index < 0 || wall_index >= (int)wall_list.size()) {
      throw RuntimeError(
          "Vertex index " + std::to_string(wall_index) + " is out of range for " + NAME_VERTEX_LIST + " of " + name + ".");
    }
  }

public:
  // extra information used in MCell4 converter
  // if this is a compartment, its parent and children are set
  // in API::set_children_compartments
  std::shared_ptr<API::GeometryObject> parent_compartment;
  GeometryObjectSet child_compartments;

  // simulation engine mapping
  partition_id_t partition_id;
  geometry_object_id_t geometry_object_id;
  vertex_index_t first_vertex_index; // index of the first vertex created in partition for this object
  wall_index_t first_wall_index;
  std::vector<vertex_index_t> vertex_indices; // vertex_list[i] has vertex index vertex_indices[i]
  std::vector<wall_index_t> wall_indices; // wall_list[i] has wall index wall_indices[i]

  BNG::compartment_id_t vol_compartment_id;
  BNG::compartment_id_t surf_compartment_id;
};

} // namespace API
} // namespace MCell

#endif // API_GEOMETRY_OBJECT_H
