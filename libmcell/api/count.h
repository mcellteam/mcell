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
  // ctor used when converting observables from BNGL
  Count(bool /*argument just to distinguish it from the generated variant*/) {
    set_all_attributes_as_default_or_unset();
    count_event = nullptr;
  }

  COUNT_CTOR()

  void postprocess_in_ctor() override {
    count_event = nullptr;
  }

  void check_semantics() const override {
    GenCount::check_semantics(); // calls also CountTerm::check_semantics
    uint num_set = get_num_set(count_expression, species_pattern, molecules_pattern, reaction_rule);
    if (num_set != 1) {
      throw ValueError(
          S("Exactly one of ") + NAME_COUNT_EXPRESSION + ", " + NAME_SPECIES_PATTERN + ", " +
          NAME_MOLECULES_PATTERN + " or " + NAME_REACTION_RULE + " must be set for " + NAME_CLASS_COUNT + ".");
    }

    if (!is_set(file_name)) {
      throw ValueError(S("Attribute ") + NAME_FILE_NAME + " must be set.");
    }

    if (is_set(count_expression)) {
      count_expression->check_that_species_or_reaction_rule_is_set();
    }

    if (every_n_timesteps < 0) {
      throw ValueError(
          S("The value of ") + NAME_EVERY_N_TIMESTEPS + " must be higher or equal to 0.");
    }
  }

  double get_current_value() override;

  std::string export_to_python(std::ostream& out, PythonExportContext& ctx) override {
    // we need to overwrite the current value for export however we do not want to change it
    // permanently
    uint64_t initial_reactions_count_orig = initial_reactions_count;
    initial_reactions_count = initial_reactions_count_export_override;
    std::string res = GenCount::export_to_python(out, ctx);
    initial_reactions_count = initial_reactions_count_orig;
    return res;
  }

  // count event, owned by Scheduler if every_n_timesteps > 0,
  // owned by World if every_n_timesteps == 0
  MolOrRxnCountEvent* count_event;
};

} // namespace API
} // namespace MCell

#endif // API_COUNT_H
