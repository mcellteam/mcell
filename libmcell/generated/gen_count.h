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

#include "api/api_common.h"
#include "api/count_term.h"


namespace MCell {
namespace API {

class Count;
class Complex;
class CountTerm;
class ReactionRule;
class Region;
class PythonExportContext;

#define COUNT_CTOR() \
    Count( \
        const std::string& name_ = STR_UNSET, \
        const std::string& file_name_ = STR_UNSET, \
        std::shared_ptr<CountTerm> count_expression_ = nullptr, \
        const double multiplier_ = 1, \
        const double every_n_timesteps_ = 1, \
        std::shared_ptr<Complex> species_pattern_ = nullptr, \
        std::shared_ptr<Complex> molecules_pattern_ = nullptr, \
        std::shared_ptr<ReactionRule> reaction_rule_ = nullptr, \
        std::shared_ptr<Region> region_ = nullptr, \
        const ExprNodeType node_type_ = ExprNodeType::LEAF, \
        std::shared_ptr<CountTerm> left_node_ = nullptr, \
        std::shared_ptr<CountTerm> right_node_ = nullptr, \
        const uint64_t initial_reactions_count_ = 0 \
    )  : GenCount(species_pattern_,molecules_pattern_,reaction_rule_,region_,node_type_,left_node_,right_node_,initial_reactions_count_) { \
      class_name = "Count"; \
      name = name_; \
      file_name = file_name_; \
      count_expression = count_expression_; \
      multiplier = multiplier_; \
      every_n_timesteps = every_n_timesteps_; \
      species_pattern = species_pattern_; \
      molecules_pattern = molecules_pattern_; \
      reaction_rule = reaction_rule_; \
      region = region_; \
      node_type = node_type_; \
      left_node = left_node_; \
      right_node = right_node_; \
      initial_reactions_count = initial_reactions_count_; \
      postprocess_in_ctor(); \
      check_semantics(); \
    } \
    Count(DefaultCtorArgType) : \
      GenCount(DefaultCtorArgType()) { \
      set_all_attributes_as_default_or_unset(); \
    }

class GenCount: public CountTerm {
public:
  GenCount( 
      std::shared_ptr<Complex> species_pattern_ = nullptr, 
      std::shared_ptr<Complex> molecules_pattern_ = nullptr, 
      std::shared_ptr<ReactionRule> reaction_rule_ = nullptr, 
      std::shared_ptr<Region> region_ = nullptr, 
      const ExprNodeType node_type_ = ExprNodeType::LEAF, 
      std::shared_ptr<CountTerm> left_node_ = nullptr, 
      std::shared_ptr<CountTerm> right_node_ = nullptr, 
      const uint64_t initial_reactions_count_ = 0 
  )  : CountTerm(species_pattern_,molecules_pattern_,reaction_rule_,region_,node_type_,left_node_,right_node_,initial_reactions_count_)  {
  }
  GenCount(DefaultCtorArgType) {
  }
  void postprocess_in_ctor() override {}
  void check_semantics() const override;
  void set_initialized() override;
  void set_all_attributes_as_default_or_unset() override;

  std::shared_ptr<Count> copy_count() const;
  std::shared_ptr<Count> deepcopy_count(py::dict = py::dict()) const;
  virtual bool __eq__(const Count& other) const;
  virtual bool eq_nonarray_attributes(const Count& other, const bool ignore_name = false) const;
  bool operator == (const Count& other) const { return __eq__(other);}
  bool operator != (const Count& other) const { return !__eq__(other);}
  std::string to_str(const bool all_details=false, const std::string ind="") const override;

  virtual std::string export_to_python(std::ostream& out, PythonExportContext& ctx);


  // --- attributes ---
  std::string file_name;
  virtual void set_file_name(const std::string& new_file_name_) {
    if (initialized) {
      throw RuntimeError("Value 'file_name' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    file_name = new_file_name_;
  }
  virtual const std::string& get_file_name() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return file_name;
  }

  std::shared_ptr<CountTerm> count_expression;
  virtual void set_count_expression(std::shared_ptr<CountTerm> new_count_expression_) {
    if (initialized) {
      throw RuntimeError("Value 'count_expression' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    count_expression = new_count_expression_;
  }
  virtual std::shared_ptr<CountTerm> get_count_expression() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return count_expression;
  }

  double multiplier;
  virtual void set_multiplier(const double new_multiplier_) {
    if (initialized) {
      throw RuntimeError("Value 'multiplier' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    multiplier = new_multiplier_;
  }
  virtual double get_multiplier() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return multiplier;
  }

  double every_n_timesteps;
  virtual void set_every_n_timesteps(const double new_every_n_timesteps_) {
    if (initialized) {
      throw RuntimeError("Value 'every_n_timesteps' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    every_n_timesteps = new_every_n_timesteps_;
  }
  virtual double get_every_n_timesteps() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return every_n_timesteps;
  }

  // --- methods ---
  virtual double get_current_value() = 0;
}; // GenCount

class Count;
py::class_<Count> define_pybinding_Count(py::module& m);
} // namespace API
} // namespace MCell

#endif // API_GEN_COUNT_H
