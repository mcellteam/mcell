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

#ifndef API_GEN_COUNT_H
#define API_GEN_COUNT_H

#include "../api/common.h"
#include "../api/count_term.h"


namespace MCell {
namespace API {

class ComplexInstance;
class CountTerm;
class ReactionRule;
class Region;

#define COUNT_CTOR() \
    Count( \
        const std::string& file_name_, \
        std::shared_ptr<CountTerm> count_expression_ = nullptr, \
        const float_t multiplier_ = 1, \
        const int every_n_timesteps_ = 1, \
        std::shared_ptr<ComplexInstance> species_pattern_ = nullptr, \
        std::shared_ptr<ComplexInstance> molecules_pattern_ = nullptr, \
        std::shared_ptr<ReactionRule> reaction_rule_ = nullptr, \
        std::shared_ptr<Region> region_ = nullptr, \
        const Orientation orientation_ = Orientation::NOT_SET, \
        const ExprNodeType node_type_ = ExprNodeType::LEAF, \
        std::shared_ptr<CountTerm> left_node_ = nullptr, \
        std::shared_ptr<CountTerm> right_node_ = nullptr \
    )  : GenCount(species_pattern_,molecules_pattern_,reaction_rule_,region_,orientation_,node_type_,left_node_,right_node_) { \
      class_name = "Count"; \
      file_name = file_name_; \
      count_expression = count_expression_; \
      multiplier = multiplier_; \
      every_n_timesteps = every_n_timesteps_; \
      species_pattern = species_pattern_; \
      molecules_pattern = molecules_pattern_; \
      reaction_rule = reaction_rule_; \
      region = region_; \
      orientation = orientation_; \
      node_type = node_type_; \
      left_node = left_node_; \
      right_node = right_node_; \
      postprocess_in_ctor();\
      check_semantics();\
    }

class GenCount: public CountTerm {
public:
  GenCount( 
      std::shared_ptr<ComplexInstance> species_pattern_ = nullptr, 
      std::shared_ptr<ComplexInstance> molecules_pattern_ = nullptr, 
      std::shared_ptr<ReactionRule> reaction_rule_ = nullptr, 
      std::shared_ptr<Region> region_ = nullptr, 
      const Orientation orientation_ = Orientation::NOT_SET, 
      const ExprNodeType node_type_ = ExprNodeType::LEAF, 
      std::shared_ptr<CountTerm> left_node_ = nullptr, 
      std::shared_ptr<CountTerm> right_node_ = nullptr 
  )  : CountTerm(species_pattern_,molecules_pattern_,reaction_rule_,region_,orientation_,node_type_,left_node_,right_node_)  {
  }
  void postprocess_in_ctor() override {}
  void check_semantics() const override;
  void set_initialized() override;
  bool __eq__(const GenCount& other) const;
  void set_all_attributes_as_default_or_unset() override;

  std::string to_str(const std::string ind="") const override;

  // --- attributes ---
  std::string file_name;
  virtual void set_file_name(const std::string& new_file_name_) {
    if (initialized) {
      throw RuntimeError("Value 'file_name' of object with name " + name + " (class " + class_name + ")"
                         "cannot be set after model was initialized.");
    }
    file_name = new_file_name_;
  }
  virtual const std::string& get_file_name() const {
    return file_name;
  }

  std::shared_ptr<CountTerm> count_expression;
  virtual void set_count_expression(std::shared_ptr<CountTerm> new_count_expression_) {
    if (initialized) {
      throw RuntimeError("Value 'count_expression' of object with name " + name + " (class " + class_name + ")"
                         "cannot be set after model was initialized.");
    }
    count_expression = new_count_expression_;
  }
  virtual std::shared_ptr<CountTerm> get_count_expression() const {
    return count_expression;
  }

  float_t multiplier;
  virtual void set_multiplier(const float_t new_multiplier_) {
    if (initialized) {
      throw RuntimeError("Value 'multiplier' of object with name " + name + " (class " + class_name + ")"
                         "cannot be set after model was initialized.");
    }
    multiplier = new_multiplier_;
  }
  virtual float_t get_multiplier() const {
    return multiplier;
  }

  int every_n_timesteps;
  virtual void set_every_n_timesteps(const int new_every_n_timesteps_) {
    if (initialized) {
      throw RuntimeError("Value 'every_n_timesteps' of object with name " + name + " (class " + class_name + ")"
                         "cannot be set after model was initialized.");
    }
    every_n_timesteps = new_every_n_timesteps_;
  }
  virtual int get_every_n_timesteps() const {
    return every_n_timesteps;
  }

  // --- methods ---
}; // GenCount

class Count;
py::class_<Count> define_pybinding_Count(py::module& m);
} // namespace API
} // namespace MCell

#endif // API_GEN_COUNT_H
