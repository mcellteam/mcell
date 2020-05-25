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

#include "../generated/gen_count.h"
#include "../api/common.h"
#include "../api/count_term.h"


namespace MCell {
namespace API {

class Count: public GenCount {
public:
  COUNT_CTOR()

  void check_semantics() const override {
    uint num_set = get_num_set(count_expression, species, reaction_rule);
    if (num_set != 1) {
      throw ValueError(
          S("Exactly one of ") + NAME_COUNT_EXPRESSION + ", " + NAME_SPECIES + " or " + NAME_REACTION_RULE +
          " must be set.");
    }

    if (is_set(count_expression)) {
      count_expression->check_that_species_or_reaction_rule_is_set();
    }
  }
};

} // namespace API
} // namespace MCell

#endif // API_COUNT_H
