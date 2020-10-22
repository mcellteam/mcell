/******************************************************************************
 *
 * Copyright (C) 2020 by
 * The Salk Institute for Biological Studies
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

#include "api/wall.h"
#include "api/geometry_object.h"

#include "world.h"
#include "partition.h"
#include "geometry.h"


using namespace std;

namespace MCell {
namespace API {

void Wall::set_is_movable(const bool new_is_movable_) {
  check_initialization();
  is_movable = new_is_movable_;

  MCell::Partition& p = world->get_partition(PARTITION_ID_INITIAL);
  MCell::Wall& w = p.get_wall(geometry_object->get_partition_wall_index(wall_index));
  w.is_movable = new_is_movable_;
}


} // namespace API
} // namespace MCell
