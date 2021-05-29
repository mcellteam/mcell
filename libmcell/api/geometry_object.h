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
#include "api/api_common.h"
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
  void postprocess_in_ctor() override;
  void check_semantics() const override;

  void translate(const std::vector<double> move) override;

  // using shorter printout when all_details is false
  std::string to_str(const bool all_details=false, const std::string ind="") const override;

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
