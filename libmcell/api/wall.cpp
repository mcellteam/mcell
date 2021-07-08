/******************************************************************************
 *
 * Copyright (C) 2020 by
 * The Salk Institute for Biological Studies
 *
 * Use of this source code is governed by an MIT-style
 * license that can be found in the LICENSE file or at
 * https://opensource.org/licenses/MIT.
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
