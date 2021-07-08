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
    res = BNG::COMPLEX_ZERO;
  }
  return res;
}


static std::string get_rxn_str(
    const std::vector<std::shared_ptr<Complex>>& reactants,
    const std::vector<std::shared_ptr<Complex>>& products,
    const double rev_rate,
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


void ReactionRule::update_reaction_rate(const BNG::rxn_rule_id_t rxn_rule_id, const double new_rate) {
  assert(is_initialized());

  if (world->scheduler.get_event_being_executed() != nullptr) {
    throw ValueError("Reaction rates cannot be changed in callbacks or in between of iterations.");
  }

  BNG::RxnRule* rxn = world->get_all_rxns().get(rxn_rule_id);
  bool updated = rxn->update_rxn_rate(new_rate);
  if (updated && rxn->is_unimol()) {
    world->reset_unimol_rxn_times(rxn_rule_id);
  }
}


void ReactionRule::set_fwd_rate(const double new_fwd_rate_) {
  if (!is_set(new_fwd_rate_)) {
    throw ValueError(S("Attribute ") + NAME_FWD_RATE + " must be set to a different value than " +
        NAME_CV_FLT_UNSET + ".");
  }

  if (is_initialized()) {
    update_reaction_rate(fwd_rxn_rule_id, new_fwd_rate_);
    cached_data_are_uptodate = false;
    fwd_rate = new_fwd_rate_;
  }
  else {
    cached_data_are_uptodate = false;
    fwd_rate = new_fwd_rate_;
  }
}


void ReactionRule::set_rev_rate(const double new_rev_rate_) {
  if (is_initialized()) {
    // TODO: add test
    if (is_set(rev_rate) && !is_set(new_rev_rate_)) {
      throw ValueError(S("Attribute ") + NAME_REV_RATE + " must be set to a different value than " +
          NAME_CV_FLT_UNSET + " is this " + NAME_CLASS_REACTION_RULE + " was initialized as a reversible reaction.");
    }

    // update the existing rule
    release_assert(rev_rxn_rule_id != BNG::RXN_RULE_ID_INVALID);
    update_reaction_rate(rev_rxn_rule_id, new_rev_rate_);
    cached_data_are_uptodate = false;
    rev_rate = new_rev_rate_;
  }
  else {
    cached_data_are_uptodate = false;
    rev_rate = new_rev_rate_;
  }
}


void ReactionRule::set_variable_rate(const std::vector<std::vector<double>> new_variable_rate_) {
  if (initialized) {
    throw RuntimeError("Value 'variable_rate' of object with name " + name + " (class " + class_name + ") "
                       "cannot be set after model was initialized.");
  }
  if (is_reversible()) {
    throw RuntimeError(S("Cannot set ") + NAME_VARIABLE_RATE + " for a reversible reaction.");
  }
  cached_data_are_uptodate = false;
  variable_rate = new_variable_rate_;
  check_variable_rate();
  // reset fwd rate so that the variable rate is used
  fwd_rate = FLT_UNSET;
}


} // namespace API
} // namespace MCell
