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

  if (!is_set(file_name)) {
    throw ValueError(S("Attribute ") + NAME_FILE_NAME + " must be set.");
  }

  if (every_n_timesteps < 0) {
    throw ValueError(
        S("The value of ") + NAME_EVERY_N_TIMESTEPS + " must be higher or equal to 0.");
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

  if (sz > gdat.size() && file_name.substr(sz - gdat_sz) == gdat) {
    output_format = CountOutputFormat::GDAT;
  }
  else if (sz > dat.size() && file_name.substr(sz - dat_sz) == dat) {
    output_format = CountOutputFormat::DAT;
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

