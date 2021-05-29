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

#ifndef API_GEN_RELEASE_PATTERN_H
#define API_GEN_RELEASE_PATTERN_H

#include "api/api_common.h"
#include "api/base_data_class.h"

namespace MCell {
namespace API {

class ReleasePattern;
class PythonExportContext;

#define RELEASE_PATTERN_CTOR() \
    ReleasePattern( \
        const std::string& name_ = STR_UNSET, \
        const double release_interval_ = TIME_INFINITY, \
        const double train_duration_ = TIME_INFINITY, \
        const double train_interval_ = TIME_INFINITY, \
        const int number_of_trains_ = 1 \
    ) { \
      class_name = "ReleasePattern"; \
      name = name_; \
      release_interval = release_interval_; \
      train_duration = train_duration_; \
      train_interval = train_interval_; \
      number_of_trains = number_of_trains_; \
      postprocess_in_ctor(); \
      check_semantics(); \
    } \
    ReleasePattern(DefaultCtorArgType) : \
      GenReleasePattern(DefaultCtorArgType()) { \
      set_all_attributes_as_default_or_unset(); \
    }

class GenReleasePattern: public BaseDataClass {
public:
  GenReleasePattern() {
  }
  GenReleasePattern(DefaultCtorArgType) {
  }
  void postprocess_in_ctor() override {}
  void check_semantics() const override;
  void set_initialized() override;
  void set_all_attributes_as_default_or_unset() override;

  std::shared_ptr<ReleasePattern> copy_release_pattern() const;
  std::shared_ptr<ReleasePattern> deepcopy_release_pattern(py::dict = py::dict()) const;
  virtual bool __eq__(const ReleasePattern& other) const;
  virtual bool eq_nonarray_attributes(const ReleasePattern& other, const bool ignore_name = false) const;
  bool operator == (const ReleasePattern& other) const { return __eq__(other);}
  bool operator != (const ReleasePattern& other) const { return !__eq__(other);}
  std::string to_str(const std::string ind="") const override;

  std::string export_to_python(std::ostream& out, PythonExportContext& ctx) override;


  // --- attributes ---
  double release_interval;
  virtual void set_release_interval(const double new_release_interval_) {
    if (initialized) {
      throw RuntimeError("Value 'release_interval' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    release_interval = new_release_interval_;
  }
  virtual double get_release_interval() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return release_interval;
  }

  double train_duration;
  virtual void set_train_duration(const double new_train_duration_) {
    if (initialized) {
      throw RuntimeError("Value 'train_duration' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    train_duration = new_train_duration_;
  }
  virtual double get_train_duration() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return train_duration;
  }

  double train_interval;
  virtual void set_train_interval(const double new_train_interval_) {
    if (initialized) {
      throw RuntimeError("Value 'train_interval' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    train_interval = new_train_interval_;
  }
  virtual double get_train_interval() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return train_interval;
  }

  int number_of_trains;
  virtual void set_number_of_trains(const int new_number_of_trains_) {
    if (initialized) {
      throw RuntimeError("Value 'number_of_trains' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    number_of_trains = new_number_of_trains_;
  }
  virtual int get_number_of_trains() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return number_of_trains;
  }

  // --- methods ---
}; // GenReleasePattern

class ReleasePattern;
py::class_<ReleasePattern> define_pybinding_ReleasePattern(py::module& m);
} // namespace API
} // namespace MCell

#endif // API_GEN_RELEASE_PATTERN_H
