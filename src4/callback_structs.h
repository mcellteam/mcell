/*
 * callback_structs.h
 *
 *  Created on: Oct 14, 2019
 *      Author: ahusar
 */

#ifndef SRC4_CALLBACK_STRUCTS_H_
#define SRC4_CALLBACK_STRUCTS_H_

#ifndef SWIG
#include "defines.h"
#endif

namespace MCell {

struct WallHitInfo {
  molecule_id_t molecule_id;
  geometry_object_id_t geometry_object_id;
  wall_id_t wall_id;
  float_t time;
  vec3_t pos;
  vec3_t pos_before_hit;
};

} /* namespace MCell */

#endif /* SRC4_CALLBACK_STRUCTS_H_ */
