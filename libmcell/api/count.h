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

#ifndef API_COUNT_H
#define API_COUNT_H

#include "generated/gen_count.h"
#include "api/api_common.h"
#include "api/count_term.h"


namespace MCell {

class MolOrRxnCountEvent;

namespace API {

class Count: public GenCount {
public:
  COUNT_CTOR()

  void postprocess_in_ctor() override {
    set_all_custom_attributes_to_default();
  }

  void set_all_custom_attributes_to_default() override {
    count_event = nullptr;
  }

  void check_semantics() const override {
    GenCount::check_semantics(); // calls also CountTerm::check_semantics
    if (!is_set(expression)) {
      throw ValueError(S("Attribute ") + NAME_EXPRESSION + " must be set.");
    }
    expression->check_that_species_or_reaction_rule_is_set();

    if (!is_set(file_name)) {
      throw ValueError(S("Attribute ") + NAME_FILE_NAME + " must be set.");
    }

    if (every_n_timesteps < 0) {
      throw ValueError(
          S("The value of ") + NAME_EVERY_N_TIMESTEPS + " must be higher or equal to 0.");
    }
  }

  double get_current_value() override;

  // count event, owned by Scheduler if every_n_timesteps > 0,
  // owned by World if every_n_timesteps == 0
  MolOrRxnCountEvent* count_event;
};

} // namespace API
} // namespace MCell

#endif // API_COUNT_H
