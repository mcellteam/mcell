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
#include "../api/base_data_class.h"
#include "../api/count_term.h"


namespace MCell {
namespace API {

class CountTerm;
class ReactionRule;
class Region;
class Species;

#define COUNT_CTOR() \
    Count( \
        const std::string& filename_, \
        std::shared_ptr<CountTerm> count_expression_ = nullptr, \
        const int every_n_timesteps_ = 1, \
        std::shared_ptr<Species> species_ = nullptr, \
        std::shared_ptr<ReactionRule> reaction_rule_ = nullptr, \
        std::shared_ptr<Region> region_ = nullptr, \
        const Orientation orientation_ = Orientation::NotSet, \
        const ExprNodeType node_type_ = ExprNodeType::Leaf, \
        std::shared_ptr<CountTerm> left_node_ = nullptr, \
        std::shared_ptr<CountTerm> right_node_ = nullptr \
    )  : GenCount(species_,reaction_rule_,region_,orientation_,node_type_,left_node_,right_node_) { \
      class_name = "Count"; \
      filename = filename_; \
      count_expression = count_expression_; \
      every_n_timesteps = every_n_timesteps_; \
      species = species_; \
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
      std::shared_ptr<Species> species_ = nullptr, 
      std::shared_ptr<ReactionRule> reaction_rule_ = nullptr, 
      std::shared_ptr<Region> region_ = nullptr, 
      const Orientation orientation_ = Orientation::NotSet, 
      const ExprNodeType node_type_ = ExprNodeType::Leaf, 
      std::shared_ptr<CountTerm> left_node_ = nullptr, 
      std::shared_ptr<CountTerm> right_node_ = nullptr 
  )  : CountTerm(species_,reaction_rule_,region_,orientation_,node_type_,left_node_,right_node_)  {
  }
  void postprocess_in_ctor() override {}
  void check_semantics() const override;
  void set_initialized() override;

  bool __eq__(const GenCount& other) const;
  std::string to_str(const std::string ind="") const override;

  // --- attributes ---
  std::string filename;
  virtual void set_filename(const std::string& new_filename_) {
    if (initialized) {
      throw RuntimeError("Value 'filename' of object with name " + name + " (class " + class_name + ")"
                         "cannot be set after model was initialized.");
    }
    filename = new_filename_;
  }
  virtual const std::string& get_filename() const {
    return filename;
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
