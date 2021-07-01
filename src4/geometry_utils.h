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

#ifndef SRC4_GEOMETRY_UTILS_H_
#define SRC4_GEOMETRY_UTILS_H_

#include "defines.h"

namespace MCell {

class Partition;
class Wall;
struct Vec3;
struct Vec2;

namespace GeometryUtils {

// some commonly used utilities for which one does not need to
// include the whole geometry_utils.inc
static inline Vec3 uv2xyz(const Vec2& a, const Wall& w, const Vec3& wall_vert0) {
  return Vec3(a.u) * w.unit_u + Vec3(a.v) * w.unit_v + wall_vert0;
}


// only the needed functions for now

static pos_t closest_interior_point(
    const Partition& p,
    const Vec3& pt,
    const Wall& w,
    Vec2& ip
);


} // namespace WallUtil

} // namespace MCell

#endif // SRC4_GEOMETRY_UTILS_H_
