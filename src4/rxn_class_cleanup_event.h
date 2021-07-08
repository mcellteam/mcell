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

#ifndef SRC4_RXN_CLASS_CLEANUP_EVENT_H_
#define SRC4_RXN_CLASS_CLEANUP_EVENT_H_

#include "base_event.h"

namespace MCell {

class World;

class RxnClassCleanupEvent: public BaseEvent {
public:
  RxnClassCleanupEvent(World* world_)
    : BaseEvent(EVENT_TYPE_INDEX_RXN_CLASS_CLEANUP),
      world(world_) {
  }

  void step() override;
  void dump(const std::string indent) const override;
private:
  World* world;
};

} /* namespace MCell */

#endif /* SRC4_RXN_CLASS_CLEANUP_EVENT_H_ */
