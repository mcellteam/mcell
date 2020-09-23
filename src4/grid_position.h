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
