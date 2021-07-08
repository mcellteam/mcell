/******************************************************************************
 *
 * Copyright (C) 2021 by
 * The Salk Institute for Biological Studies
 *
 * Use of this source code is governed by an MIT-style
 * license that can be found in the LICENSE file or at
 * https://opensource.org/licenses/MIT.
 *
******************************************************************************/

#ifndef API_GEN_COUNT_H
#define API_GEN_COUNT_H

#include "api/api_common.h"
#include "api/base_data_class.h"

namespace MCell {
namespace API {

class Count;
class CountTerm;
class PythonExportContext;

#define COUNT_CTOR() \
    Count( \
        const std::string& name_ = STR_UNSET, \
        const std::string& file_name_ = STR_UNSET, \
        std::shared_ptr<CountTerm> expression_ = nullptr, \
        const double multiplier_ = 1, \
        const double every_n_timesteps_ = 1, \
        const CountOutputFormat output_format_ = CountOutputFormat::AUTOMATIC_FROM_EXTENSION \
    ) { \
      class_name = "Count"; \
      name = name_; \
      file_name = file_name_; \
      expression = expression_; \
      multiplier = multiplier_; \
      every_n_timesteps = every_n_timesteps_; \
      output_format = output_format_; \
      postprocess_in_ctor(); \
      check_semantics(); \
    } \
    Count(DefaultCtorArgType) : \
      GenCount(DefaultCtorArgType()) { \
      set_all_attributes_as_default_or_unset(); \
      set_all_custom_attributes_to_default(); \
    }

class GenCount: public BaseDataClass {
public:
  GenCount() {
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

  std::string export_to_python(std::ostream& out, PythonExportContext& ctx) override;


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

  std::shared_ptr<CountTerm> expression;
  virtual void set_expression(std::shared_ptr<CountTerm> new_expression_) {
    if (initialized) {
      throw RuntimeError("Value 'expression' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    expression = new_expression_;
  }
  virtual std::shared_ptr<CountTerm> get_expression() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return expression;
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

  CountOutputFormat output_format;
  virtual void set_output_format(const CountOutputFormat new_output_format_) {
    if (initialized) {
      throw RuntimeError("Value 'output_format' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    output_format = new_output_format_;
  }
  virtual CountOutputFormat get_output_format() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return output_format;
  }

  // --- methods ---
  virtual double get_current_value() = 0;
}; // GenCount

class Count;
py::class_<Count> define_pybinding_Count(py::module& m);
} // namespace API
} // namespace MCell

#endif // API_GEN_COUNT_H
