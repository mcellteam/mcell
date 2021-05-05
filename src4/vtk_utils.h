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
#include "geometry.h"

namespace MCell {

class World;
class Partition;
class GeometryObject;

namespace VtkUtils {


class GeomObjectInfo {
public:
  GeomObjectInfo(const partition_id_t partition_id_, const geometry_object_id_t geometry_object_id_)
    : for_counted_objects(true), partition_id(partition_id_), geometry_object_id(geometry_object_id_) {
  }

  GeomObjectInfo(const std::string& name_)
    : for_counted_objects(false), partition_id(PARTITION_ID_INVALID), geometry_object_id(GEOMETRY_OBJECT_ID_INVALID),
      name(name_) {
  }

  // true - for counted objects
  // false - for compartments
  bool for_counted_objects;

  // partition_id and geometry_object_id are used when
  // we are computing containment for counted objects
  // TODO: remove partition_id - it is not needed because we can get this info from geometry_object_id
  partition_id_t partition_id;
  geometry_object_id_t geometry_object_id;

  // name is used when we are computing containment for compartments
  std::string name;

  vtkSmartPointer<vtkPolyData> polydata;

  GeometryObject& get_geometry_object_noconst(World* world) const;
  const GeometryObject& get_geometry_object(const World* world) const;

  // comparison uses geometry_object_id only, used in maps
  bool operator < (const GeomObjectInfo& other) const {
    if (for_counted_objects) {
      return geometry_object_id < other.geometry_object_id;
    }
    else {
      return name < other.name;
    }
  }

  bool operator == (const GeomObjectInfo& other) const {
    if (for_counted_objects) {
      return geometry_object_id == other.geometry_object_id;
    }
    else {
      return name == other.name;
    }
  }
};


typedef std::vector<GeomObjectInfo> GeomObjectInfoVector;

// containment mapping of counted geometry objects
typedef std::map<GeomObjectInfo, std::set<GeomObjectInfo>> ContainmentMap;
typedef std::set<GeomObjectInfo> IntersectingSet;


// return true if counted volumes were correctly set up
bool initialize_counted_volumes(World* world, bool& has_intersecting_counted_objects);

#if 0
bool is_point_inside_counted_volume(GeometryObject& obj, const Vec3& point);
#endif

// world is nullptr for compartment hierarchy computation
bool compute_containement_mapping(
    const World* world, const GeomObjectInfoVector& counted_objects,
    ContainmentMap& contained_in_mapping,
    IntersectingSet& intersecting_objects);

const GeomObjectInfo* get_direct_parent_info(
    const GeomObjectInfo& obj_info, const ContainmentMap& contained_in_mapping);


// auxiliary function to compute volume, not related to counted volumes but uses VTK
double get_geometry_object_volume(const World* world, const GeometryObject& obj);

void export_geometry_objects_to_obj(
    const World* world, const GeometryObjectVector& objs, const std::string& file_prefix);

}; // namespace VtkUtils

} // namespace MCell

#endif // SRC4_COUNTED_VOLUME_UTIL_H_
