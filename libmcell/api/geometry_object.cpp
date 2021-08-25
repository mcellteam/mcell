/******************************************************************************
 *
 * Copyright (C) 2021 by
 * The Salk Institute for Biological Studies
 *
 * Use of this source code is governed by an MIT-style
 * license that can be found in the LICENSE file or at
 * https://opensource.org/licenses/MIT.
 *
******************************************************************************/

#include "api/geometry_object.h"

#include "api/surface_class.h"

using namespace std;

namespace MCell {
namespace API {

void GeometryObject::postprocess_in_ctor() {
  set_all_custom_attributes_to_default();
}


void GeometryObject::set_all_custom_attributes_to_default() {
  Region::set_all_custom_attributes_to_default();

  // overwrite value in Region construction
  region_type = RegionType::VOLUME;

  node_type = RegionNodeType::LEAF_GEOMETRY_OBJECT;
  partition_id = PARTITION_ID_INVALID;
  geometry_object_id = GEOMETRY_OBJECT_ID_INVALID;
  first_vertex_index = VERTEX_INDEX_INVALID;
  parent_compartment = nullptr;
  vol_compartment_id = BNG::COMPARTMENT_ID_INVALID;
  surf_compartment_id = BNG::COMPARTMENT_ID_INVALID;
}


void GeometryObject::check_semantics() const {
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


void GeometryObject::translate(const std::vector<double> move) {
  if (move.size() != 3) {
    throw ValueError(S("Argument ") + NAME_MOVE + " must be a list containing exactly 3 floats.");
  }
  if (geometry_object_id != GEOMETRY_OBJECT_ID_INVALID) {
    throw RuntimeError(S("Method ") + NAME_TRANSLATE + " may be called only before model initialization.");
  }
  for (auto& v: vertex_list) {
    v[0] += move[0];
    v[1] += move[1];
    v[2] += move[2];
  }
}


std::string GeometryObject::to_str(const bool all_details, const std::string ind) const {
  if (!all_details) {
    std::stringstream ss;
    ss << get_object_name() << ": " <<
        "name=" << name << ", " <<
        "is_bngl_compartment=" << is_bngl_compartment << ", " <<
        "surface_compartment_name=" << surface_compartment_name << ", " <<
        "surface_class=" << "(" << ((surface_class != nullptr) ? surface_class->to_str(all_details, ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  ";

    return ss.str();
  }
  else {
    return GenGeometryObject::to_str(true, ind);
  }
}


} // namespace API
} // namespace MCell
