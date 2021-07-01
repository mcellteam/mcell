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

#ifndef API_COUNT_TERM_H
#define API_COUNT_TERM_H

#include "generated/gen_count_term.h"
#include "api/api_common.h"
#include "api/region.h"
#include "api/reaction_rule.h"
#include "api/complex.h"

namespace MCell {
namespace API {

class CountTerm: public GenCountTerm, public std::enable_shared_from_this<CountTerm> {
public:
  COUNT_TERM_CTOR()

  void postprocess_in_ctor() override {
    set_all_custom_attributes_to_default();
  }

  void set_all_custom_attributes_to_default() override {
    initial_reactions_count_export_override = 0;
  }

  void check_semantics() const override {
    GenCountTerm::check_semantics();

    if (is_set(reaction_rule) && is_set(reaction_rule->rev_rate)) {
      throw ValueError(S("Reversible reactions cannot counted because it is not clear which direction should be counted.") +
          " Error for " + NAME_CLASS_COUNT_TERM + " " + name + ". Split the reaction rule into its forward and reverse variants if needed.");
    }

    // TODO: add test for this case, not sure if the check works correctly because pattern does not have to be initialized
    if (is_set(get_pattern(false)) &&
        BNG::get_in_or_out_compartment_id(get_pattern()->compartment_name) != BNG::COMPARTMENT_ID_INVALID) {
      throw ValueError(
          S(NAME_CLASS_COUNT) + " or " + NAME_CLASS_COUNT_TERM + " must not use compartment class name " +
          get_pattern()->get_primary_compartment_name() + ".");
    }
  }

  // called from Count::check_semantics()
  void check_that_species_or_reaction_rule_is_set() {

    if (node_type == ExprNodeType::LEAF) {
      uint num_set = get_num_set(species_pattern, molecules_pattern, reaction_rule);
      if (num_set != 1) {
        // NOTE: does not give much information on where to search for the error
        throw ValueError(
            S("Exactly one of ") + NAME_SPECIES_PATTERN + ", " +
            NAME_MOLECULES_PATTERN + " or " + NAME_REACTION_RULE + " must be set for one of the " +
            NAME_CLASS_COUNT_TERM + " used in " + NAME_CLASS_COUNT + ".");
      }
    }
    else {
      left_node->check_that_species_or_reaction_rule_is_set();
      right_node->check_that_species_or_reaction_rule_is_set();
    }

  }

  std::shared_ptr<CountTerm> create_expr_term(ExprNodeType op, std::shared_ptr<CountTerm> op2) {
    std::shared_ptr<CountTerm> res = std::make_shared<CountTerm>(DefaultCtorArgType());
    res->node_type = op;
    res->left_node = shared_from_this();
    res->right_node = op2;
    return res;
  }

  std::shared_ptr<CountTerm> __add__(std::shared_ptr<CountTerm> op2) override {
    return create_expr_term(ExprNodeType::ADD, op2);
  }

  std::shared_ptr<CountTerm> __sub__(std::shared_ptr<CountTerm> op2) override {
    return create_expr_term(ExprNodeType::SUB, op2);
  }

  // manually added, may return empty shared_ptr if reaction_rule
  // is counted
  std::shared_ptr<Complex> get_pattern(bool must_be_present = true) const {
    if (is_set(species_pattern)) {
      return species_pattern;
    }
    else if (is_set(molecules_pattern)) {
      return molecules_pattern;
    }
    else {
      assert(!must_be_present || is_set(reaction_rule));
      return std::shared_ptr<Complex>(nullptr);
    }
  }


  std::string export_to_python(std::ostream& out, PythonExportContext& ctx) override {
    // we need to overwrite the current value for export however we do not want to change it
    // permanently
    uint64_t initial_reactions_count_orig = initial_reactions_count;
    initial_reactions_count = initial_reactions_count_export_override;
    std::string res = GenCountTerm::export_to_python(out, ctx);
    initial_reactions_count = initial_reactions_count_orig;
    return res;
  }

  uint64_t initial_reactions_count_export_override;
};

} // namespace API
} // namespace MCell

#endif // API_COUNT_TERM_H
