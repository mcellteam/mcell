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
