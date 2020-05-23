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

#ifndef API_GEN_COUNT_TERM_H
#define API_GEN_COUNT_TERM_H

#include "../api/common.h"
#include "../api/base_data_class.h"

namespace MCell {
namespace API {

class CountTerm;
class ReactionRule;
class Region;
class Species;

#define COUNT_TERM_CTOR() \
    CountTerm( \
        std::shared_ptr<Species> species_ = nullptr, \
        std::shared_ptr<ReactionRule> reaction_rule_ = nullptr, \
        std::shared_ptr<Region> region_ = nullptr, \
        const ExprNodeType node_type_ = ExprNodeType::Leaf, \
        std::shared_ptr<CountTerm> left_node_ = nullptr, \
        std::shared_ptr<CountTerm> right_node_ = nullptr \
    ) { \
      class_name = "CountTerm"; \
      species = species_; \
      reaction_rule = reaction_rule_; \
      region = region_; \
      node_type = node_type_; \
      left_node = left_node_; \
      right_node = right_node_; \
      postprocess_in_ctor();\
      check_semantics();\
    }

class GenCountTerm: public BaseDataClass {
public:
  void postprocess_in_ctor() override {}
  void check_semantics() const override;
  void set_initialized() override;

  bool __eq__(const GenCountTerm& other) const;
  std::string to_str(const std::string ind="") const override;

  // --- attributes ---
  std::shared_ptr<Species> species;
  virtual void set_species(std::shared_ptr<Species> new_species_) {
    if (initialized) {
      throw RuntimeError("Value 'species' of object with name " + name + " (class " + class_name + ")"
                         "cannot be set after model was initialized.");
    }
    species = new_species_;
  }
  virtual std::shared_ptr<Species> get_species() const {
    return species;
  }

  std::shared_ptr<ReactionRule> reaction_rule;
  virtual void set_reaction_rule(std::shared_ptr<ReactionRule> new_reaction_rule_) {
    if (initialized) {
      throw RuntimeError("Value 'reaction_rule' of object with name " + name + " (class " + class_name + ")"
                         "cannot be set after model was initialized.");
    }
    reaction_rule = new_reaction_rule_;
  }
  virtual std::shared_ptr<ReactionRule> get_reaction_rule() const {
    return reaction_rule;
  }

  std::shared_ptr<Region> region;
  virtual void set_region(std::shared_ptr<Region> new_region_) {
    if (initialized) {
      throw RuntimeError("Value 'region' of object with name " + name + " (class " + class_name + ")"
                         "cannot be set after model was initialized.");
    }
    region = new_region_;
  }
  virtual std::shared_ptr<Region> get_region() const {
    return region;
  }

  ExprNodeType node_type;
  virtual void set_node_type(const ExprNodeType new_node_type_) {
    if (initialized) {
      throw RuntimeError("Value 'node_type' of object with name " + name + " (class " + class_name + ")"
                         "cannot be set after model was initialized.");
    }
    node_type = new_node_type_;
  }
  virtual ExprNodeType get_node_type() const {
    return node_type;
  }

  std::shared_ptr<CountTerm> left_node;
  virtual void set_left_node(std::shared_ptr<CountTerm> new_left_node_) {
    if (initialized) {
      throw RuntimeError("Value 'left_node' of object with name " + name + " (class " + class_name + ")"
                         "cannot be set after model was initialized.");
    }
    left_node = new_left_node_;
  }
  virtual std::shared_ptr<CountTerm> get_left_node() const {
    return left_node;
  }

  std::shared_ptr<CountTerm> right_node;
  virtual void set_right_node(std::shared_ptr<CountTerm> new_right_node_) {
    if (initialized) {
      throw RuntimeError("Value 'right_node' of object with name " + name + " (class " + class_name + ")"
                         "cannot be set after model was initialized.");
    }
    right_node = new_right_node_;
  }
  virtual std::shared_ptr<CountTerm> get_right_node() const {
    return right_node;
  }

  // --- methods ---
  virtual std::shared_ptr<CountTerm> __add__(std::shared_ptr<CountTerm> op2) = 0;
  virtual std::shared_ptr<CountTerm> __sub__(std::shared_ptr<CountTerm> op2) = 0;
}; // GenCountTerm

class CountTerm;
py::class_<CountTerm> define_pybinding_CountTerm(py::module& m);
} // namespace API
} // namespace MCell

#endif // API_GEN_COUNT_TERM_H
