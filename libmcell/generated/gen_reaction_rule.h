/******************************************************************************
 *
 * Copyright (C) 2020 by
 * The Salk Institute for Biological Studies and
 * Pittsburgh Supercomputing Center, Carnegie Mellon University
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

#ifndef API_GEN_REACTION_RULE_H
#define API_GEN_REACTION_RULE_H

#include "../api/common.h"

namespace MCell {
namespace API {

class ComplexInstance;

#define REACTION_RULE_CTOR() \
    ReactionRule( \
        const std::string& name_ = STR_UNSET, \
        const std::vector<std::shared_ptr<ComplexInstance>> reactants_ = std::vector<std::shared_ptr<ComplexInstance>>(), \
        const std::vector<std::shared_ptr<ComplexInstance>> products_ = std::vector<std::shared_ptr<ComplexInstance>>(), \
        const float_t fwd_rate_ = FLT_UNSET, \
        const float_t rev_rate_ = FLT_UNSET \
    ) { \
      class_name = "ReactionRule"; \
      name = name_; \
      reactants = reactants_; \
      products = products_; \
      fwd_rate = fwd_rate_; \
      rev_rate = rev_rate_; \
      postprocess_in_ctor();\
      check_semantics();\
    }

class GenReactionRule: public BaseDataClass {
public:
  void postprocess_in_ctor() override {}
  void check_semantics() const override;
  std::string to_str(const std::string ind="") const override;

  // --- attributes ---
  std::vector<std::shared_ptr<ComplexInstance>> reactants;
  virtual void set_reactants(const std::vector<std::shared_ptr<ComplexInstance>> new_reactants_) {
    reactants = new_reactants_;
  }
  virtual std::vector<std::shared_ptr<ComplexInstance>> get_reactants() const {
    return reactants;
  }

  std::vector<std::shared_ptr<ComplexInstance>> products;
  virtual void set_products(const std::vector<std::shared_ptr<ComplexInstance>> new_products_) {
    products = new_products_;
  }
  virtual std::vector<std::shared_ptr<ComplexInstance>> get_products() const {
    return products;
  }

  float_t fwd_rate;
  virtual void set_fwd_rate(const float_t new_fwd_rate_) {
    fwd_rate = new_fwd_rate_;
  }
  virtual float_t get_fwd_rate() const {
    return fwd_rate;
  }

  float_t rev_rate;
  virtual void set_rev_rate(const float_t new_rev_rate_) {
    rev_rate = new_rev_rate_;
  }
  virtual float_t get_rev_rate() const {
    return rev_rate;
  }

  // --- methods ---
}; // GenReactionRule

class ReactionRule;
py::class_<ReactionRule> define_pybinding_ReactionRule(py::module& m);
} // namespace API
} // namespace MCell

#endif // API_GEN_REACTION_RULE_H
