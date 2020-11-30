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


#ifndef SRC4_SPECIES_FLAGS_ANALYZER_H_
#define SRC4_SPECIES_FLAGS_ANALYZER_H_

#include "bng/base_flag.h"

namespace BNG {
class Species;
}

namespace MCell {

class BaseEvent;
class MolOrRxnCountEvent;

class SpeciesFlagsAnalyzer: public BNG::BaseCustomFlagsAnalyzer {
public:
  SpeciesFlagsAnalyzer()
    : initialized(false) {
  }

  // vector base_count_events may be empty
  void initialize(
      std::vector<BaseEvent*>& scheduled_count_events,
      std::vector<MolOrRxnCountEvent*> unscheduled_count_events);

  // returns a mask of all custom flags for species,
  // used to clear any pre-existing flag values
  uint get_custom_species_flags_mask() const override {
    assert(initialized);
    return BNG::SPECIES_FLAG_NEEDS_COUNTED_VOLUME;
  }

  // returns a mask of all custom flags for species that should be set
  // currently figures out whether the species need SPECIES_FLAG_NEEDS_COUNTED_VOLUME
  uint get_custom_species_flags_to_set(const BNG::Species& species) const override;

private:
  bool initialized;
  std::vector<MolOrRxnCountEvent*> count_events;
};

} /* namespace MCell */

#endif /* SRC4_SPECIES_FLAGS_ANALYZER_H_ */
