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

#ifndef API_GEN_REACTION_RULE_H
#define API_GEN_REACTION_RULE_H

#include "../api/common.h"
#include "../api/base_data_class.h"

namespace MCell {
namespace API {

class ComplexInstance;

#define REACTION_RULE_CTOR() \
    ReactionRule( \
        const std::string& name_ = STR_UNSET, \
        const std::vector<std::shared_ptr<ComplexInstance>> reactants_ = std::vector<std::shared_ptr<ComplexInstance>>(), \
        const std::vector<std::shared_ptr<ComplexInstance>> products_ = std::vector<std::shared_ptr<ComplexInstance>>(), \
        const float_t fwd_rate_ = FLT_UNSET, \
        const std::string& rev_name_ = STR_UNSET, \
        const float_t rev_rate_ = FLT_UNSET, \
        const std::vector<std::vector<float_t>> variable_rate_ = std::vector<std::vector<float_t>>() \
    ) { \
      class_name = "ReactionRule"; \
      name = name_; \
      reactants = reactants_; \
      products = products_; \
      fwd_rate = fwd_rate_; \
      rev_name = rev_name_; \
      rev_rate = rev_rate_; \
      variable_rate = variable_rate_; \
      postprocess_in_ctor();\
      check_semantics();\
    }

class GenReactionRule: public BaseDataClass {
public:
  void postprocess_in_ctor() override {}
  void check_semantics() const override;
  void set_initialized() override;
  bool __eq__(const GenReactionRule& other) const;
  void set_all_attributes_as_default_or_unset() override;

  std::string to_str(const std::string ind="") const override;

  // --- attributes ---
  std::vector<std::shared_ptr<ComplexInstance>> reactants;
  virtual void set_reactants(const std::vector<std::shared_ptr<ComplexInstance>> new_reactants_) {
    if (initialized) {
      throw RuntimeError("Value 'reactants' of object with name " + name + " (class " + class_name + ")"
                         "cannot be set after model was initialized.");
    }
    reactants = new_reactants_;
  }
  virtual std::vector<std::shared_ptr<ComplexInstance>> get_reactants() const {
    return reactants;
  }

  std::vector<std::shared_ptr<ComplexInstance>> products;
  virtual void set_products(const std::vector<std::shared_ptr<ComplexInstance>> new_products_) {
    if (initialized) {
      throw RuntimeError("Value 'products' of object with name " + name + " (class " + class_name + ")"
                         "cannot be set after model was initialized.");
    }
    products = new_products_;
  }
  virtual std::vector<std::shared_ptr<ComplexInstance>> get_products() const {
    return products;
  }

  float_t fwd_rate;
  virtual void set_fwd_rate(const float_t new_fwd_rate_) {
    if (initialized) {
      throw RuntimeError("Value 'fwd_rate' of object with name " + name + " (class " + class_name + ")"
                         "cannot be set after model was initialized.");
    }
    fwd_rate = new_fwd_rate_;
  }
  virtual float_t get_fwd_rate() const {
    return fwd_rate;
  }

  std::string rev_name;
  virtual void set_rev_name(const std::string& new_rev_name_) {
    if (initialized) {
      throw RuntimeError("Value 'rev_name' of object with name " + name + " (class " + class_name + ")"
                         "cannot be set after model was initialized.");
    }
    rev_name = new_rev_name_;
  }
  virtual const std::string& get_rev_name() const {
    return rev_name;
  }

  float_t rev_rate;
  virtual void set_rev_rate(const float_t new_rev_rate_) {
    if (initialized) {
      throw RuntimeError("Value 'rev_rate' of object with name " + name + " (class " + class_name + ")"
                         "cannot be set after model was initialized.");
    }
    rev_rate = new_rev_rate_;
  }
  virtual float_t get_rev_rate() const {
    return rev_rate;
  }

  std::vector<std::vector<float_t>> variable_rate;
  virtual void set_variable_rate(const std::vector<std::vector<float_t>> new_variable_rate_) {
    if (initialized) {
      throw RuntimeError("Value 'variable_rate' of object with name " + name + " (class " + class_name + ")"
                         "cannot be set after model was initialized.");
    }
    variable_rate = new_variable_rate_;
  }
  virtual std::vector<std::vector<float_t>> get_variable_rate() const {
    return variable_rate;
  }

  // --- methods ---
}; // GenReactionRule

class ReactionRule;
py::class_<ReactionRule> define_pybinding_ReactionRule(py::module& m);
} // namespace API
} // namespace MCell

#endif // API_GEN_REACTION_RULE_H
