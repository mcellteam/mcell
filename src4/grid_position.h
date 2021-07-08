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

#ifndef SRC4_GRID_POSITION_H_
#define SRC4_GRID_POSITION_H_

#include "defines.h"

namespace MCell {

class Partition;
class GridPos;

namespace GridPosition {

Vec2 find_closest_position(const Partition& p, const GridPos& gp1, const GridPos& gp2);

} // namespace GridPosition
} // namespace MCell


#endif // SRC4_GRID_POSITION_H_
