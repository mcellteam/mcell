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

namespace MCell {

class World;
class GeometryObject;

namespace CountedVolumesUtil {

// return true if counted volumes were correctly set up
bool initialize_counted_volumes(World* world, bool& has_intersecting_counted_objects);

bool is_point_inside_counted_volume(GeometryObject& obj, const Vec3& point);

};

} // namespace MCell

#endif // SRC4_COUNTED_VOLUME_UTIL_H_
