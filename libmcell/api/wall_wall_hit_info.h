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

#ifndef API_WALL_WALL_HIT_INFO_H
#define API_WALL_WALL_HIT_INFO_H

#include "generated/gen_wall_wall_hit_info.h"
#include "api/api_common.h"

namespace MCell {
namespace API {

class WallWallHitInfo: public GenWallWallHitInfo {
public:
  WALL_WALL_HIT_INFO_CTOR_NOARGS()
};

} // namespace API
} // namespace MCell

#endif // API_WALL_WALL_HIT_INFO_H
