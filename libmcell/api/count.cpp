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

#include "api/count.h"

#include "mol_or_rxn_count_event.h"

namespace MCell {
namespace API {

float_t Count::get_current_value() {
  if (count_event == nullptr) {
    throw RuntimeError(S(NAME_CLASS_COUNT) + " with name " + name + " was not initialized.");
  }

  if (is_set(reaction_rule)) {
    throw RuntimeError(S("Calling of ") + NAME_GET_CURRENT_VALUE + " to count reactions is not supported " +
        "( counted reaction rule: " + reaction_rule->to_bngl_str() + ").");
  }

  return count_event->get_single_count_value(count_event_index);
}

} // namespace API
} // namespace MCell

