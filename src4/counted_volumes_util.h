/******************************************************************************
 *
 * Copyright (C) 2020 by
 * The Salk Institute for Biological Studies and
 * Pittsburgh Supercomputing Center, Carnegie Mellon University
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

#ifndef SRC4_COUNTED_VOLUME_UTIL_H_
#define SRC4_COUNTED_VOLUME_UTIL_H_

#include <vector>
#include <map>

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>

#include "defines.h"

namespace MCell {

class World;
class GeometryObject;

namespace CountedVolumesUtil {


struct GeomObjectInfo {

  GeomObjectInfo(const partition_id_t partition_id_, const geometry_object_id_t geometry_object_id_)
    : partition_id(partition_id_), geometry_object_id(geometry_object_id_) {
  }

  // we are using IDs because in the future we might need to create new
  // geometry objects and pointers would become invalidated
  partition_id_t partition_id;
  geometry_object_id_t geometry_object_id;

  vtkSmartPointer<vtkPolyData> polydata;

  GeometryObject& get_geometry_object_noconst(World* world) const;
  const GeometryObject& get_geometry_object(const World* world) const;

  // comparison uses geometry_object_id only, used in maps
  bool operator < (const GeomObjectInfo& second) const {
    return geometry_object_id < second.geometry_object_id;
  }

  bool operator == (const GeomObjectInfo& second) const {
    return geometry_object_id == second.geometry_object_id;
  }
};


typedef std::vector<GeomObjectInfo> GeomObjectInfoVector;

// containment mapping of counted geometry objects
typedef std::map<GeomObjectInfo, uint_set<GeomObjectInfo>> ContainmentMap;
typedef uint_set<geometry_object_id_t> IntersectingSet;


// return true if counted volumes were correctly set up
bool initialize_counted_volumes(World* world, bool& has_intersecting_counted_objects);

bool is_point_inside_counted_volume(GeometryObject& obj, const Vec3& point);

bool compute_containement_mapping(
    const GeomObjectInfoVector& counted_objects,
    ContainmentMap& contained_in_mapping,
    IntersectingSet& intersecting_objects);

};

} // namespace MCell

#endif // SRC4_COUNTED_VOLUME_UTIL_H_
