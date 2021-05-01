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

#include "api/api_common.h"
#include "api/base_data_class.h"

namespace MCell {
namespace API {

class ReactionRule;
class Complex;
class PythonExportContext;

#define REACTION_RULE_CTOR() \
    ReactionRule( \
        const std::string& name_ = STR_UNSET, \
        const std::vector<std::shared_ptr<Complex>> reactants_ = std::vector<std::shared_ptr<Complex>>(), \
        const std::vector<std::shared_ptr<Complex>> products_ = std::vector<std::shared_ptr<Complex>>(), \
        const double fwd_rate_ = FLT_UNSET, \
        const std::string& rev_name_ = STR_UNSET, \
        const double rev_rate_ = FLT_UNSET, \
        const std::vector<std::vector<double>> variable_rate_ = std::vector<std::vector<double>>(), \
        const bool is_intermembrane_surface_reaction_ = false \
    ) { \
      class_name = "ReactionRule"; \
      name = name_; \
      reactants = reactants_; \
      products = products_; \
      fwd_rate = fwd_rate_; \
      rev_name = rev_name_; \
      rev_rate = rev_rate_; \
      variable_rate = variable_rate_; \
      is_intermembrane_surface_reaction = is_intermembrane_surface_reaction_; \
      postprocess_in_ctor();\
      check_semantics();\
    }

class GenReactionRule: public BaseDataClass {
public:
  void postprocess_in_ctor() override {}
  void check_semantics() const override;
  void set_initialized() override;
  void set_all_attributes_as_default_or_unset() override;

  virtual bool __eq__(const ReactionRule& other) const;
  virtual bool eq_nonarray_attributes(const ReactionRule& other, const bool ignore_name = false) const;
  bool operator == (const ReactionRule& other) const { return __eq__(other);}
  bool operator != (const ReactionRule& other) const { return !__eq__(other);}
  std::string to_str(const std::string ind="") const override;

  std::string export_to_python(std::ostream& out, PythonExportContext& ctx) override;
  virtual std::string export_vec_reactants(std::ostream& out, PythonExportContext& ctx, const std::string& parent_name);
  virtual std::string export_vec_products(std::ostream& out, PythonExportContext& ctx, const std::string& parent_name);
  virtual std::string export_vec_variable_rate(std::ostream& out, PythonExportContext& ctx, const std::string& parent_name);


  // --- attributes ---
  std::vector<std::shared_ptr<Complex>> reactants;
  virtual void set_reactants(const std::vector<std::shared_ptr<Complex>> new_reactants_) {
    if (initialized) {
      throw RuntimeError("Value 'reactants' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    reactants = new_reactants_;
  }
  virtual std::vector<std::shared_ptr<Complex>>& get_reactants() {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return reactants;
  }

  std::vector<std::shared_ptr<Complex>> products;
  virtual void set_products(const std::vector<std::shared_ptr<Complex>> new_products_) {
    if (initialized) {
      throw RuntimeError("Value 'products' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    products = new_products_;
  }
  virtual std::vector<std::shared_ptr<Complex>>& get_products() {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return products;
  }

  double fwd_rate;
  virtual void set_fwd_rate(const double new_fwd_rate_) {
    if (initialized) {
      throw RuntimeError("Value 'fwd_rate' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    fwd_rate = new_fwd_rate_;
  }
  virtual double get_fwd_rate() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return fwd_rate;
  }

  std::string rev_name;
  virtual void set_rev_name(const std::string& new_rev_name_) {
    if (initialized) {
      throw RuntimeError("Value 'rev_name' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    rev_name = new_rev_name_;
  }
  virtual const std::string& get_rev_name() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return rev_name;
  }

  double rev_rate;
  virtual void set_rev_rate(const double new_rev_rate_) {
    if (initialized) {
      throw RuntimeError("Value 'rev_rate' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    rev_rate = new_rev_rate_;
  }
  virtual double get_rev_rate() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return rev_rate;
  }

  std::vector<std::vector<double>> variable_rate;
  virtual void set_variable_rate(const std::vector<std::vector<double>> new_variable_rate_) {
    if (initialized) {
      throw RuntimeError("Value 'variable_rate' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    variable_rate = new_variable_rate_;
  }
  virtual std::vector<std::vector<double>>& get_variable_rate() {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return variable_rate;
  }

  bool is_intermembrane_surface_reaction;
  virtual void set_is_intermembrane_surface_reaction(const bool new_is_intermembrane_surface_reaction_) {
    if (initialized) {
      throw RuntimeError("Value 'is_intermembrane_surface_reaction' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    is_intermembrane_surface_reaction = new_is_intermembrane_surface_reaction_;
  }
  virtual bool get_is_intermembrane_surface_reaction() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return is_intermembrane_surface_reaction;
  }

  // --- methods ---
  virtual std::string to_bngl_str() const = 0;
}; // GenReactionRule

class ReactionRule;
py::class_<ReactionRule> define_pybinding_ReactionRule(py::module& m);
} // namespace API
} // namespace MCell

#endif // API_GEN_REACTION_RULE_H
