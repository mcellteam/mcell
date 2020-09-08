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

#ifndef SRC4_SPECIES_CLEANUP_EVENT_H_
#define SRC4_SPECIES_CLEANUP_EVENT_H_

#include "base_event.h"

namespace MCell {

class World;

/**
 * Removes all surface and volume species that whose instantiation count is 0.
 * Removes all rxn classes because they might reference the removed species.
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
  World* world;
};

} /* namespace MCell */

#endif /* SRC4_SPECIES_CLEANUP_EVENT_H_ */
