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


#ifndef SRC4_REACTION_H_
#define SRC4_REACTION_H_

#include "defines.h"

namespace mcell {

typedef int32_t orientation_t;
const orientation_t ORIENTATION_DOWN = -1;
const orientation_t ORIENTATION_NONE = 0;
const orientation_t ORIENTATION_UP = 1;

struct species_with_orientation_t {
  species_with_orientation_t()
    : species_id(SPECIES_ID_INVALID), orientation(ORIENTATION_NONE) {
  }
  species_with_orientation_t(
      const species_id_t species_id_,
      const orientation_t orientation_)
    : species_id(species_id_), orientation(orientation_) {
  }

  species_id_t species_id;
  orientation_t orientation;
};


class reaction_t {
public:
	// TODO: add comments
  std::string name;
  float_t rate_constant;
  float_t max_fixed_p;
  float_t min_noreaction_p;
  std::vector<species_with_orientation_t> reactants;
  std::vector<species_with_orientation_t> products;

  void dump(const std::string ind) const;
};


/**
 * Used as a pair molecule id, remaining timestep for molecules newly created in diffusion.
 * Using name action instead of event because events are handled by scheduler and are ordered by time.
 * These actions are simply processes in a queue (FIFO) manner.
 *
 * Used in diffuse_react _event_t and in partition_t.
 */
class diffuse_or_unimol_react_action_t {
public:
  enum type_t {
    DIFFUSE,
    UNIMOL_REACT
  };

  diffuse_or_unimol_react_action_t(
      const molecule_id_t id_,
      const float_t scheduled_time_,
      const type_t type_,
      const reaction_t* unimol_rx_ = nullptr)
    :
      id(id_),
      scheduled_time(scheduled_time_),
      type(type_),
      unimol_rx(unimol_rx_) {
    if (type == UNIMOL_REACT) {
      assert(unimol_rx_ != nullptr);
    }
  }

  // defined because of usage in calendar_t
  const diffuse_or_unimol_react_action_t& operator->() const {
     return *this;
  }

  molecule_id_t id;
  float_t scheduled_time; // this is the scheduled time
  type_t type;
  const reaction_t* unimol_rx; // when type is UNIMOL_REACT
};

} // namespace mcell

#endif // SRC4_REACTION_H_
