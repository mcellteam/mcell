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

#include "api/complex.h"
#include "api/reaction_rule.h"
#include "bng/bngl_names.h"

#include "world.h"

#include <set>

using namespace std;

namespace MCell {
namespace API {

bool ReactionRule::__eq__(const ReactionRule& other) const {
  // ignore name when comparing
  if (!eq_nonarray_attributes(other, true)) {
    return false;
  }

  if (variable_rate != other.variable_rate) {
    return false;
  }

  string n1 = get_canonical_name();
  string n2 = other.get_canonical_name();

  return n1 == n2;
}


static std::string get_rxn_side_str(
    const std::vector<std::shared_ptr<Complex>>& cplxs,
    const bool canonical) {

  string res;
  if (!cplxs.empty()) {
    for (size_t i = 0; i < cplxs.size(); i++) {
      if (!canonical) {
        res += cplxs[i]->to_bngl_str();
      }
      else {
        res += cplxs[i]->get_canonical_name();
      }
      if (i + 1 != cplxs.size()) {
        res += " + ";
      }
    }
  }
  else {
    res = BNG::COMPLEX_Null;
  }
  return res;
}


static std::string get_rxn_str(
    const std::vector<std::shared_ptr<Complex>>& reactants,
    const std::vector<std::shared_ptr<Complex>>& products,
    const float_t rev_rate,
    bool canonical) {

  string res;
  res = get_rxn_side_str(reactants, canonical);

  if (is_set(rev_rate)) {
    res += " <-> ";
  }
  else {
    res += " -> ";
  }

  res += get_rxn_side_str(products, canonical);
  return res;
}


std::string ReactionRule::to_bngl_str_w_orientation(
    bool replace_orientation_w_up_down_compartments) const {

  return get_rxn_str(reactants, products, rev_rate, replace_orientation_w_up_down_compartments);
}


static void sort_by_canonical_name(std::vector<std::shared_ptr<Complex>>& vec) {
  std::sort(vec.begin(), vec.end(),
      [](const std::shared_ptr<Complex>& a, const std::shared_ptr<Complex>& b) -> bool {
          return a->get_canonical_name() < b->get_canonical_name();
      });
}


std::string ReactionRule::get_canonical_name() const {

  std::vector<std::shared_ptr<Complex>> r_sorted = reactants;
  sort_by_canonical_name(r_sorted);

  std::vector<std::shared_ptr<Complex>> p_sorted = products;
  sort_by_canonical_name(p_sorted);

  return get_rxn_str(r_sorted, p_sorted, rev_rate, true);
}


void ReactionRule::set_fwd_rate(const float_t new_fwd_rate_) {
  if (!is_set(new_fwd_rate_)) {
    throw ValueError(S("Attribute ") + NAME_FWD_RATE + " must be set to a different value than " +
        NAME_CV_FLT_UNSET + ".");
  }

  if (is_initialized()) {
    // update the existing rule
    BNG::RxnRule* rxn = world->get_all_rxns().get(fwd_rxn_rule_id);
    rxn->update_rxn_rate(new_fwd_rate_);
    cached_data_are_uptodate = false;
    fwd_rate = new_fwd_rate_;
  }
  else {
    cached_data_are_uptodate = false;
    fwd_rate = new_fwd_rate_;
  }
}

void ReactionRule::set_rev_rate(const float_t new_rev_rate_) {
  if (is_initialized()) {
    // TODO: add test
    if (is_set(rev_rate) && !is_set(new_rev_rate_)) {
      throw ValueError(S("Attribute ") + NAME_REV_RATE + " must be set to a different value than " +
          NAME_CV_FLT_UNSET + " is this " + NAME_CLASS_REACTION_RULE + " was initialized as a reversible reaction.");
    }

    // update the existing rule
    release_assert(rev_rxn_rule_id != BNG::RXN_RULE_ID_INVALID);
    BNG::RxnRule* rxn = world->get_all_rxns().get(rev_rxn_rule_id);
    rxn->update_rxn_rate(new_rev_rate_);
    cached_data_are_uptodate = false;
    rev_rate = new_rev_rate_;
  }
  else {
    cached_data_are_uptodate = false;
    rev_rate = new_rev_rate_;
  }
}

} // namespace API
} // namespace MCell
