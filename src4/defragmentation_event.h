/******************************************************************************
 *
 * Copyright (C) 2019 by
 * The Salk Institute for Biological Studies and
 * Pittsburgh Supercomputing Center, Carnegie Mellon University
 *
 * Use of this source code is governed by an MIT-style
 * license that can be found in the LICENSE file or at
 * https://opensource.org/licenses/MIT.
 *
******************************************************************************/

#ifndef SRC4_DEFRAGMENTATION_EVENT_H_
#define SRC4_DEFRAGMENTATION_EVENT_H_

#include "base_event.h"

namespace MCell {

/**
 * When a reaction occurs, a hole is created in the parition's molecules array,
 * defragmentation "squeezes" this array so that there are no defuct molecules anymore;
 * updates all related data
 */
class DefragmentationEvent: public BaseEvent {
public:
  DefragmentationEvent(World* world_)
    : BaseEvent(EVENT_TYPE_INDEX_DEFRAGMENTATION),
      world(world_) {
  }

  void step() override;
  void dump(const std::string indent) const override;
private:
  World* world;
};

} // namespace mcell

#endif // SRC4_DEFRAGMENTATION_EVENT_H_
