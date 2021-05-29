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

#include "api/api_common.h"
#include "api/base_data_class.h"

namespace MCell {
namespace API {

class Complex;
class CountTerm;
class ReactionRule;
class Region;
class PythonExportContext;

#define COUNT_TERM_CTOR() \
    CountTerm( \
        std::shared_ptr<Complex> species_pattern_ = nullptr, \
        std::shared_ptr<Complex> molecules_pattern_ = nullptr, \
        std::shared_ptr<ReactionRule> reaction_rule_ = nullptr, \
        std::shared_ptr<Region> region_ = nullptr, \
        const ExprNodeType node_type_ = ExprNodeType::LEAF, \
        std::shared_ptr<CountTerm> left_node_ = nullptr, \
        std::shared_ptr<CountTerm> right_node_ = nullptr, \
        const uint64_t initial_reactions_count_ = 0 \
    ) { \
      class_name = "CountTerm"; \
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
    CountTerm(DefaultCtorArgType) : \
      GenCountTerm(DefaultCtorArgType()) { \
      set_all_attributes_as_default_or_unset(); \
    }

class GenCountTerm: public BaseDataClass {
public:
  GenCountTerm() {
  }
  GenCountTerm(DefaultCtorArgType) {
  }
  void postprocess_in_ctor() override {}
  void check_semantics() const override;
  void set_initialized() override;
  void set_all_attributes_as_default_or_unset() override;

  std::shared_ptr<CountTerm> copy_count_term() const;
  std::shared_ptr<CountTerm> deepcopy_count_term(py::dict = py::dict()) const;
  virtual bool __eq__(const CountTerm& other) const;
  virtual bool eq_nonarray_attributes(const CountTerm& other, const bool ignore_name = false) const;
  bool operator == (const CountTerm& other) const { return __eq__(other);}
  bool operator != (const CountTerm& other) const { return !__eq__(other);}
  std::string to_str(const std::string ind="") const override;

  std::string export_to_python(std::ostream& out, PythonExportContext& ctx) override;


  // --- attributes ---
  std::shared_ptr<Complex> species_pattern;
  virtual void set_species_pattern(std::shared_ptr<Complex> new_species_pattern_) {
    if (initialized) {
      throw RuntimeError("Value 'species_pattern' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    species_pattern = new_species_pattern_;
  }
  virtual std::shared_ptr<Complex> get_species_pattern() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return species_pattern;
  }

  std::shared_ptr<Complex> molecules_pattern;
  virtual void set_molecules_pattern(std::shared_ptr<Complex> new_molecules_pattern_) {
    if (initialized) {
      throw RuntimeError("Value 'molecules_pattern' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    molecules_pattern = new_molecules_pattern_;
  }
  virtual std::shared_ptr<Complex> get_molecules_pattern() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return molecules_pattern;
  }

  std::shared_ptr<ReactionRule> reaction_rule;
  virtual void set_reaction_rule(std::shared_ptr<ReactionRule> new_reaction_rule_) {
    if (initialized) {
      throw RuntimeError("Value 'reaction_rule' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    reaction_rule = new_reaction_rule_;
  }
  virtual std::shared_ptr<ReactionRule> get_reaction_rule() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return reaction_rule;
  }

  std::shared_ptr<Region> region;
  virtual void set_region(std::shared_ptr<Region> new_region_) {
    if (initialized) {
      throw RuntimeError("Value 'region' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    region = new_region_;
  }
  virtual std::shared_ptr<Region> get_region() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return region;
  }

  ExprNodeType node_type;
  virtual void set_node_type(const ExprNodeType new_node_type_) {
    if (initialized) {
      throw RuntimeError("Value 'node_type' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    node_type = new_node_type_;
  }
  virtual ExprNodeType get_node_type() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return node_type;
  }

  std::shared_ptr<CountTerm> left_node;
  virtual void set_left_node(std::shared_ptr<CountTerm> new_left_node_) {
    if (initialized) {
      throw RuntimeError("Value 'left_node' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    left_node = new_left_node_;
  }
  virtual std::shared_ptr<CountTerm> get_left_node() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return left_node;
  }

  std::shared_ptr<CountTerm> right_node;
  virtual void set_right_node(std::shared_ptr<CountTerm> new_right_node_) {
    if (initialized) {
      throw RuntimeError("Value 'right_node' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    right_node = new_right_node_;
  }
  virtual std::shared_ptr<CountTerm> get_right_node() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return right_node;
  }

  uint64_t initial_reactions_count;
  virtual void set_initial_reactions_count(const uint64_t new_initial_reactions_count_) {
    if (initialized) {
      throw RuntimeError("Value 'initial_reactions_count' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    initial_reactions_count = new_initial_reactions_count_;
  }
  virtual uint64_t get_initial_reactions_count() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return initial_reactions_count;
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
