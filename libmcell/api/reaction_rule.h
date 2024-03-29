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

#ifndef API_REACTION_RULE_H
#define API_REACTION_RULE_H

#include "generated/gen_reaction_rule.h"
#include "api/api_common.h"
#include "bng/bng.h"

namespace MCell {

class World;

namespace API {

class ReactionRule: public GenReactionRule {
public:
  REACTION_RULE_CTOR()

  void postprocess_in_ctor() override {
    set_all_custom_attributes_to_default();
  }

  void set_all_custom_attributes_to_default() override {
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
        throw ValueError(S("Variable rates cannot be set along with fwd_rate or rev_rate in the constructor of ") +
            NAME_CLASS_REACTION_RULE + ".");
      }

      check_variable_rate();
    }
  }

  void check_variable_rate() const {
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

  bool __eq__(const ReactionRule& other) const override;

  std::string to_bngl_str() const override {
    return to_bngl_str_w_orientation();
  }

  void set_fwd_rate(const double new_fwd_rate_) override;
  void set_rev_rate(const double new_rev_rate_) override;
  void set_variable_rate(const std::vector<std::vector<double>> new_variable_rate_) override;

  // added methods
  bool is_reversible() const {
    return is_set(rev_rate);
  }

  bool eq_reactants_and_products(const ReactionRule& other) const;

  std::string to_bngl_str_w_orientation(bool replace_orientation_w_up_down_compartments = false) const;

  std::string get_canonical_name() const;

  bool warn_if_adding_identical_object() const { return true; }

  // simulation engine mapping
  bool is_initialized() const {
    assert((fwd_rxn_rule_id == BNG::RXN_RULE_ID_INVALID) == (world == nullptr) &&
        "When initialized, both values must be set");
    return world != nullptr;
  }
  BNG::rxn_rule_id_t fwd_rxn_rule_id;
  BNG::rxn_rule_id_t rev_rxn_rule_id;
  World* world;

private:
  void update_reaction_rate(const BNG::rxn_rule_id_t rxn_rule_id, const double new_rate);
};

} // namespace API
} // namespace MCell

#endif // API_REACTION_RULE_H
