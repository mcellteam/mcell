
/******************************************************************************
 *
 * Copyright (C) 2021 by
 * The Salk Institute for Biological Studies and
 * Pittsburgh Supercomputing Center, Carnegie Mellon University
 *
 * Use of this source code is governed by an MIT-style
 * license that can be found in the LICENSE file or at
 * https://opensource.org/licenses/MIT.
 *
******************************************************************************/

#ifndef SRC4_MOL_ORDER_SHUFFLE_EVENT_H_
#define SRC4_MOL_ORDER_SHUFFLE_EVENT_H_

#include "base_event.h"

namespace MCell {

class MolOrderShuffleEvent: public BaseEvent {
public:
  MolOrderShuffleEvent(World* world_)
    : BaseEvent(EVENT_TYPE_INDEX_MOL_SHUFFLE),
      world(world_) {
  }

  void step() override;
  void dump(const std::string indent) const override;
private:
  World* world;
};

} // namespace mcell

#endif // SRC4_MOL_ORDER_SHUFFLE_EVENT_H_
