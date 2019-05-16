/******************************************************************************
 *
 * Copyright (C) 2006-2017 by
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

#ifndef SRC4_GEOMETRY_UTILS_INC_
#define SRC4_GEOMETRY_UTILS_INC_

/**
 * This file is directly included into diffuse_react_event.cpp.
 * The reason why this is not a standard .cpp + .h file is to gove the compiler
 * the opportunity to inline these functions into methods of diffuse&react event.
 */
#include <vector>

#include "diffuse_react_event.h"
#include "defines.h"
#include "world.h"
#include "partition.h"
#include "geometry.h"
#include "debug_config.h"

namespace mcell {

namespace geom_util {

/***************************************************************************
wall_bounding_box:
  In: a wall
      vector to store one corner of the bounding box for that wall
      vector to store the opposite corner
  Out: No return value.  The vectors are set to define the smallest box
       that contains the wall.
***************************************************************************/
static inline void get_wall_bounding_box(
    const vec3_t w_vert[VERTICES_IN_TRIANGLE],
    vec3_t& llf, vec3_t& urb
) {
  llf.x = urb.x = w_vert[0].x;
  llf.y = urb.y = w_vert[0].y;
  llf.z = urb.z = w_vert[0].z;

  if (w_vert[1].x < llf.x)
    llf.x = w_vert[1].x;
  else if (w_vert[1].x > urb.x)
    urb.x = w_vert[1].x;
  if (w_vert[2].x < llf.x)
    llf.x = w_vert[2].x;
  else if (w_vert[2].x > urb.x)
    urb.x = w_vert[2].x;

  if (w_vert[1].y < llf.y)
    llf.y = w_vert[1].y;
  else if (w_vert[1].y > urb.y)
    urb.y = w_vert[1].y;
  if (w_vert[2].y < llf.y)
    llf.y = w_vert[2].y;
  else if (w_vert[2].y > urb.y)
    urb.y = w_vert[2].y;

  if (w_vert[1].z < llf.z)
    llf.z = w_vert[1].z;
  else if (w_vert[1].z > urb.z)
    urb.z = w_vert[1].z;
  if (w_vert[2].z < llf.z)
    llf.z = w_vert[2].z;
  else if (w_vert[2].z > urb.z)
    urb.z = w_vert[2].z;
}

} // namespace geom_util

} // namespace mcell

#endif // SRC4_GEOMETRY_UTILS_INC_
