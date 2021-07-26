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

#include "api/count.h"

#include "mol_or_rxn_count_event.h"

namespace MCell {
namespace API {

void Count::postprocess_in_ctor() {
  set_all_custom_attributes_to_default();

  // determine count type
  set_automatically_output_format_if_needed();
}


void Count::check_semantics() const {
  GenCount::check_semantics(); // calls also CountTerm::check_semantics
  if (!is_set(expression)) {
    throw ValueError(S("Attribute ") + NAME_EXPRESSION + " must be set.");
  }
  expression->check_that_species_or_reaction_rule_is_set();

  if (every_n_timesteps < 0) {
    throw ValueError(
        S("The value of ") + NAME_EVERY_N_TIMESTEPS + " must be higher or equal to 0.");
  }

  if (output_format == CountOutputFormat::UNSET) {
    throw ValueError(S("Attribute ") + NAME_OUTPUT_FORMAT + " must be set.");
  }
}


void Count::set_automatically_output_format_if_needed() {

  if (output_format != CountOutputFormat::AUTOMATIC_FROM_EXTENSION) {
    // specific value set
    return;
  }

  const std::string gdat = ".gdat";
  size_t gdat_sz = gdat.size();
  const std::string dat = ".dat";
  size_t dat_sz = dat.size();
  size_t sz = file_name.size();

  if (!is_set(file_name) || (sz > dat.size() && file_name.substr(sz - dat_sz) == dat)) {
    output_format = CountOutputFormat::DAT;
  }
  else if (sz > gdat.size() && file_name.substr(sz - gdat_sz) == gdat) {
    output_format = CountOutputFormat::GDAT;
  }
  else {
    throw ValueError(S("Cannot automatically determine ") + NAME_OUTPUT_FORMAT + ", " + NAME_FILE_NAME +
        " must have either .dat or .gdat extension for automatic detection.");
  }
}


double Count::get_current_value() {
  if (count_event == nullptr) {
    throw RuntimeError(S(NAME_CLASS_COUNT) + " with name " + name + " was not initialized.");
  }

  return count_event->get_single_count_value();
}

} // namespace API
} // namespace MCell

