/*
 * callback_structs.h
 *
 *  Created on: Oct 14, 2019
 *      Author: ahusar
 */

#ifndef SRC4_LEGACY_CALLBACK_STRUCTS_H_
#define SRC4_LEGACY_CALLBACK_STRUCTS_H_

#ifndef SWIG
#include "defines.h"
#endif

namespace MCell {

struct LegacyWallHitInfo {
  molecule_id_t molecule_id;
  geometry_object_id_t geometry_object_id;
  wall_id_t wall_id;
  float_t time;
  Vec3 pos;
  float_t time_before_hit;
  Vec3 pos_before_hit;
};

} /* namespace MCell */

#endif /* SRC4_LEGACY_CALLBACK_STRUCTS_H_ */
