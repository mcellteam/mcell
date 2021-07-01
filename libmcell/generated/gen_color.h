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

#ifndef API_GEN_COLOR_H
#define API_GEN_COLOR_H

#include "api/api_common.h"
#include "api/base_data_class.h"

namespace MCell {
namespace API {

class Color;
class PythonExportContext;

#define COLOR_CTOR() \
    Color( \
        const double red_ = FLT_UNSET, \
        const double green_ = FLT_UNSET, \
        const double blue_ = FLT_UNSET, \
        const double alpha_ = 1, \
        const uint rgba_ = 0 \
    ) { \
      class_name = "Color"; \
      red = red_; \
      green = green_; \
      blue = blue_; \
      alpha = alpha_; \
      rgba = rgba_; \
      postprocess_in_ctor(); \
      check_semantics(); \
    } \
    Color(DefaultCtorArgType) : \
      GenColor(DefaultCtorArgType()) { \
      set_all_attributes_as_default_or_unset(); \
      set_all_custom_attributes_to_default(); \
    }

class GenColor: public BaseDataClass {
public:
  GenColor() {
  }
  GenColor(DefaultCtorArgType) {
  }
  void postprocess_in_ctor() override {}
  void check_semantics() const override;
  void set_initialized() override;
  void set_all_attributes_as_default_or_unset() override;

  std::shared_ptr<Color> copy_color() const;
  std::shared_ptr<Color> deepcopy_color(py::dict = py::dict()) const;
  virtual bool __eq__(const Color& other) const;
  virtual bool eq_nonarray_attributes(const Color& other, const bool ignore_name = false) const;
  bool operator == (const Color& other) const { return __eq__(other);}
  bool operator != (const Color& other) const { return !__eq__(other);}
  std::string to_str(const bool all_details=false, const std::string ind="") const override;

  std::string export_to_python(std::ostream& out, PythonExportContext& ctx) override;


  // --- attributes ---
  double red;
  virtual void set_red(const double new_red_) {
    if (initialized) {
      throw RuntimeError("Value 'red' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    red = new_red_;
  }
  virtual double get_red() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return red;
  }

  double green;
  virtual void set_green(const double new_green_) {
    if (initialized) {
      throw RuntimeError("Value 'green' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    green = new_green_;
  }
  virtual double get_green() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return green;
  }

  double blue;
  virtual void set_blue(const double new_blue_) {
    if (initialized) {
      throw RuntimeError("Value 'blue' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    blue = new_blue_;
  }
  virtual double get_blue() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return blue;
  }

  double alpha;
  virtual void set_alpha(const double new_alpha_) {
    if (initialized) {
      throw RuntimeError("Value 'alpha' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    alpha = new_alpha_;
  }
  virtual double get_alpha() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return alpha;
  }

  uint rgba;
  virtual void set_rgba(const uint new_rgba_) {
    if (initialized) {
      throw RuntimeError("Value 'rgba' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    rgba = new_rgba_;
  }
  virtual uint get_rgba() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return rgba;
  }

  // --- methods ---
}; // GenColor

class Color;
py::class_<Color> define_pybinding_Color(py::module& m);
} // namespace API
} // namespace MCell

#endif // API_GEN_COLOR_H
