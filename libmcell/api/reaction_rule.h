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

#ifndef API_REACTION_RULE_H
#define API_REACTION_RULE_H

#include "generated/gen_reaction_rule.h"
#include "api/common.h"
#include "bng/bng.h"

namespace MCell {
namespace API {

class ReactionRule: public GenReactionRule {
public:
  REACTION_RULE_CTOR()

  void postprocess_in_ctor() override {
    fwd_rxn_rule_id = BNG::RXN_RULE_ID_INVALID;
    rev_rxn_rule_id = BNG::RXN_RULE_ID_INVALID;
  }

  void check_semantics() const override {
    GenReactionRule::check_semantics();

    if (is_set(name) && is_set(rev_rate)) {
      if (!is_set(rev_name)) {
        throw ValueError(
            S("If name of a reaction is set, reversible reaction must have its ") + NAME_REV_NAME + " set as well."
            " Error for " + name + "."
        );
      }
    }

    if (is_set(rev_name) && !is_set(rev_rate)) {
      throw ValueError(
          S("Parameter ") + NAME_REV_NAME + " must not be set when " + NAME_REV_RATE + " is not set."
          " Error for " + name + "."
      );
    }

    if (is_set(variable_rate)) {
      if (is_set(fwd_rate) || is_set(rev_rate)) {
        throw ValueError("Variable rates cannot be set along with fwd_rate or rev_rate.");
      }

      for (auto& time_and_rate: variable_rate) {
        if (time_and_rate.size() != 2) {
          std::string msg = "[]";
          if (time_and_rate.size() >= 1) {
            msg = " item with time " + std::to_string(time_and_rate[0]);
          }
          throw ValueError("Variable rate array must contain pairs [time, rate], error for " + msg + ".");
        }
      }
    }
  }

  bool __eq__(const ReactionRule& other) const override;

  std::string to_bngl_str() const override;

  // added methods
  bool eq_reactants_and_products(const ReactionRule& other) const;


  // FIXME:
  std::string get_canonical_name() const { return to_bngl_str(); }

  // simulation engine mapping
  BNG::rxn_rule_id_t fwd_rxn_rule_id;
  BNG::rxn_rule_id_t rev_rxn_rule_id;
};

} // namespace API
} // namespace MCell

#endif // API_REACTION_RULE_H
