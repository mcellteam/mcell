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

#ifndef SRC4_COUNTED_VOLUME_UTIL_H_
#define SRC4_COUNTED_VOLUME_UTIL_H_

#include <vector>
#include <map>
#include <set>
#include <memory>
#include <string>

#include "defines.h"
#include "geometry.h"
#include "geometry_analysis.hpp"

namespace MCell {

class World;
class Partition;
class GeometryObject;

namespace CountedVolumeUtils {


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

  // the object's cleaned closed mesh, prepared for containment tests
  std::shared_ptr<const GeometryAnalysis::Solid> solid;

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

// lengths of the longest containment paths, filled by get_direct_parent_info;
// one cache serves all calls for the same containment mapping
typedef std::map<GeomObjectInfo, int> PathLengthCache;


// return true if counted volumes were correctly set up
bool initialize_counted_volumes(World* world, bool& has_intersecting_counted_objects);

// world is nullptr for compartment hierarchy computation
bool compute_containement_mapping(
    const World* world, const GeomObjectInfoVector& counted_objects,
    ContainmentMap& contained_in_mapping,
    IntersectingSet& intersecting_objects);

const GeomObjectInfo* get_direct_parent_info(
    const GeomObjectInfo& obj_info, const ContainmentMap& contained_in_mapping,
    PathLengthCache& path_lengths);


// auxiliary function to compute volume, not related to counted volumes,
// returns FLT_INVALID if the object is not closed
double get_geometry_object_volume(const World* world, const GeometryObject& obj);

void export_geometry_objects_to_obj(
    const World* world, const GeometryObjectVector& objs, const std::string& file_prefix);

}; // namespace CountedVolumeUtils

} // namespace MCell

#endif // SRC4_COUNTED_VOLUME_UTIL_H_
