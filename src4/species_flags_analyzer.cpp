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

#include "species_flags_analyzer.h"

#include "defines.h"
#include "mol_or_rxn_count_event.h"

namespace MCell {

void SpeciesFlagsAnalyzer::initialize(std::vector<BaseEvent*>& base_count_events) {

  for (BaseEvent* event: base_count_events) {
    assert(event->type_index == EVENT_TYPE_INDEX_MOL_OR_RXN_COUNT);
    MolOrRxnCountEvent* count_event = dynamic_cast<MolOrRxnCountEvent*>(event);
    assert(count_event != nullptr);
    count_events.push_back(count_event);
  }
  initialized = true;
}


uint SpeciesFlagsAnalyzer::get_custom_species_flags_to_set(const BNG::Species& species) const {
  bool needs_counted_volume = false;
  for (MolOrRxnCountEvent* count_event: count_events) {
    if (count_event->species_needs_counted_volume(species.id)) {
      needs_counted_volume = true;
      break;
    }
  }
  if (needs_counted_volume) {
    return BNG::SPECIES_FLAG_NEEDS_COUNTED_VOLUME;
  }
  else {
    return 0; // no flags to be set
  }
}


} /* namespace MCell */
