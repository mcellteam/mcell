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

bool are_coplanar(const Partition& p, const Wall& w1, const Wall& w2, const pos_t eps);
bool are_coincident(const Partition& p, const Wall& w1, const Wall& w2, const pos_t eps);

// w1 and w2 are assumed to be coplanar
bool coplanar_walls_overlap(const Partition& p, const Wall& w1, const Wall& w2);

} // namespace WallOverlap
} // namespace MCell

#endif // SRC4_WALL_OVERLAP_H_
