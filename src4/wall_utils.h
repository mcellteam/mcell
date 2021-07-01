/******************************************************************************
 *
 * Copyright (C) 2020 by
 * The Salk Institute for Biological Studies and
 * Pittsburgh Supercomputing Center, Carnegie Mellon University
 *
 * Use of this source code is governed by an MIT-style
 * license that can be found in the LICENSE file or at
 * https://opensource.org/licenses/MIT.
 *
******************************************************************************/

#ifndef SRC4_WALL_UTILS_H_
#define SRC4_WALL_UTILS_H_

namespace MCell {

class Partition;
class Wall;
struct Vec3;

namespace WallUtils {

// only the needed functions for now

static int wall_in_box(
    const Partition& p,
    const Wall& w,
    const Vec3& llf, const Vec3& urb
);


} // namespace WallUtil

} // namespace MCell

#endif // SRC4_WALL_UTILS_H_

