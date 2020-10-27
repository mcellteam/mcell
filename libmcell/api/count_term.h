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

#ifndef API_COUNT_TERM_H
#define API_COUNT_TERM_H

#include "generated/gen_count_term.h"
#include "api/common.h"
#include "api/region.h"
#include "api/reaction_rule.h"
#include "api/complex.h"

namespace MCell {
namespace API {

class CountTerm: public GenCountTerm, public std::enable_shared_from_this<CountTerm> {
public:
  COUNT_TERM_CTOR()

  void check_semantics() const override {
    GenCountTerm::check_semantics();
    if (is_set(region)) {
      if (region->node_type != RegionNodeType::LEAF_GEOMETRY_OBJECT && region->node_type != RegionNodeType::LEAF_SURFACE_REGION) {
        throw ValueError(S("Only simple regions of type ") + NAME_CLASS_GEOMETRY_OBJECT + " and " + NAME_CLASS_REGION +
            " are supported now. Error for " + NAME_CLASS_COUNT_TERM + " " + name);
      }
    }

    if (is_set(reaction_rule) && is_set(reaction_rule->rev_rate)) {
      throw ValueError(S("Reversible reactions cannot counted because it is not clear which direction should be counted.") +
          " Error for " + NAME_CLASS_COUNT_TERM + " " + name + ". Split the reaction rule into its forward and reverse variants if needed.");
    }


    if (is_set(get_pattern()) &&
        BNG::get_in_or_out_compartment_id(get_pattern()->compartment_name) != BNG::COMPARTMENT_ID_INVALID) {
      throw ValueError(
          S(NAME_CLASS_COUNT) + " or " + NAME_CLASS_COUNT_TERM + " must not use compartment class name " +
          get_pattern()->compartment_name + ".");
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
    std::shared_ptr<CountTerm> res = std::make_shared<CountTerm>();
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
  std::shared_ptr<Complex> get_pattern() const {
    if (is_set(species_pattern)) {
      return species_pattern;
    }
    else if (is_set(molecules_pattern)) {
      return molecules_pattern;
    }
    else {
      assert(is_set(reaction_rule));
      return std::shared_ptr<Complex>(nullptr);
    }
  }
};

} // namespace API
} // namespace MCell

#endif // API_COUNT_TERM_H
