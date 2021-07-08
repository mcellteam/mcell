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

#ifndef SRC4_SPECIES_CLEANUP_EVENT_H_
#define SRC4_SPECIES_CLEANUP_EVENT_H_

#include "base_event.h"

namespace MCell {

class World;

/**
 * Removes all surface and volume species that whose instantiation count is 0.
 * Removes all rxn classes because they might reference the removed species.
 * Also removes unused reactant classes.
 */
class SpeciesCleanupEvent: public BaseEvent {
public:
  SpeciesCleanupEvent(World* world_)
    : BaseEvent(EVENT_TYPE_INDEX_RXN_CLASS_CLEANUP),
      world(world_) {
  }

  void step() override;
  void dump(const std::string indent) const override;
private:

  void remove_unused_reactant_classes();

  World* world;
};

} /* namespace MCell */

#endif /* SRC4_SPECIES_CLEANUP_EVENT_H_ */
