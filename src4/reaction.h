/******************************************************************************
 *
 * Copyright (C) 2019 by
 * The Salk Institute for Biological Studies and
 * Pittsburgh Supercomputing Center, Carnegie Mellon University
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

// TODO: rename to rxn.h/cpp

#ifndef SRC4_REACTION_H_
#define SRC4_REACTION_H_

#include "bng/bng.h"

#include "defines.h"

namespace MCell {

class SimulationConfig;
class Partition;


/**
 * Used as a pair molecule id, remaining timestep for molecules newly created in diffusion.
 * Using name action instead of event because events are handled by scheduler and are ordered by time.
 * These actions are simply processes in a queue (FIFO) manner.
 *
 * Used in diffuse_react _event_t and in partition_t.
 */
class DiffuseOrUnimolRxnAction {
public:
  enum class Type {
    DIFFUSE,
    UNIMOL_REACT
  };

  // DIFFUSE action
  DiffuseOrUnimolRxnAction(
      const DiffuseOrUnimolRxnAction::Type type_,
      const molecule_id_t id_,
      const float_t scheduled_time_,
      const WallTileIndexPair& where_created_this_iteration_)
    :
      id(id_),
      scheduled_time(scheduled_time_),
      type(type_),
      unimol_rx(nullptr),
      where_created_this_iteration(where_created_this_iteration_) {

    assert(scheduled_time >= 0.0);
    assert(type == Type::DIFFUSE);
    // position where the molecule was created may be invalid when it was not a result of surface reaction
  }

  // UNIMOL_REACT action
  DiffuseOrUnimolRxnAction(
      const DiffuseOrUnimolRxnAction::Type type_,
      const molecule_id_t id_,
      const float_t scheduled_time_,
      const BNG::RxnClass* unimol_rx_)
    :
      id(id_),
      scheduled_time(scheduled_time_),
      type(type_),
      unimol_rx(unimol_rx_) {
    assert(scheduled_time >= 0.0);
    assert(type == Type::UNIMOL_REACT);
    assert(unimol_rx != nullptr);
  }

  // defined because of usage in calendar_t
  const DiffuseOrUnimolRxnAction& operator->() const {
     return *this;
  }

  molecule_id_t id;
  float_t scheduled_time; // this is the scheduled time
  Type type;

  // when type is UNIMOL_REACT
  const BNG::RxnClass* unimol_rx;

  // when type is DIFFUSE
  // used to avoid rebinding for surf+vol->surf+vol reactions
  WallTileIndexPair where_created_this_iteration;
};

} // namespace mcell

#endif // SRC4_REACTION_H_
