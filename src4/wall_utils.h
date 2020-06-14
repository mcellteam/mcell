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

#ifndef SRC4_WALL_UTILS_H_
#define SRC4_WALL_UTILS_H_

namespace MCell {

class Partition;
class Wall;
class Vec3;

namespace WallUtil {

// only the needed functions for now

static int wall_in_box(
    const Partition& p,
    const Wall& w,
    const Vec3& llf, const Vec3& urb
);


} // namespace WallUtil

} // namespace MCell

#endif // SRC4_WALL_UTILS_H_

