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

#ifndef API_WALL_H
#define API_WALL_H

#include "generated/gen_wall.h"
#include "api/api_common.h"

namespace MCell {
namespace API {

class Wall: public GenWall {
public:
  WALL_CTOR_NOARGS()

  // -- overrides ---
  void set_is_movable(const bool new_is_movable_) override;
};

} // namespace API
} // namespace MCell

#endif // API_WALL_H
