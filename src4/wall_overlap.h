/******************************************************************************
 *
 * Copyright (C) 2006-2017,2020 by
 * The Salk Institute for Biological Studies and
 * Pittsburgh Supercomputing Center, Carnegie Mellon University
 *
 * Use of this source code is governed by an MIT-style
 * license that can be found in the LICENSE file or at
 * https://opensource.org/licenses/MIT.
 *
******************************************************************************/

#ifndef SRC4_WALL_OVERLAP_H_
#define SRC4_WALL_OVERLAP_H_

#include "defines.h"

namespace MCell {
class Wall;

namespace WallOverlap {

// rand_vec is a value of three random values used when sorting walls for overlap detection
// (TODO: not completely sure why a random value is needed there)
bool check_for_overlapped_walls(Partition& p, const Vec3& rand_vec);

} // namespace WallOverlap
} // namespace MCell

#endif // SRC4_WALL_OVERLAP_H_
