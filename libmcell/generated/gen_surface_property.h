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

#ifndef API_GEN_SURFACE_PROPERTY_H
#define API_GEN_SURFACE_PROPERTY_H

#include "api/common.h"
#include "api/base_data_class.h"

namespace MCell {
namespace API {

class SurfaceProperty;
class Complex;

#define SURFACE_PROPERTY_CTOR() \
    SurfaceProperty( \
        const SurfacePropertyType type_ = SurfacePropertyType::UNSET, \
        std::shared_ptr<Complex> affected_complex_pattern_ = nullptr, \
        const float_t concentration_ = FLT_UNSET \
    ) { \
      class_name = "SurfaceProperty"; \
      type = type_; \
      affected_complex_pattern = affected_complex_pattern_; \
      concentration = concentration_; \
      postprocess_in_ctor();\
      check_semantics();\
    }

class GenSurfaceProperty: public BaseDataClass {
public:
  void postprocess_in_ctor() override {}
  void check_semantics() const override;
  void set_initialized() override;
  void set_all_attributes_as_default_or_unset() override;

  virtual bool __eq__(const SurfaceProperty& other) const;
  bool operator == (const SurfaceProperty& other) const { return __eq__(other);}
  std::string to_str(const std::string ind="") const override;

  // --- attributes ---
  SurfacePropertyType type;
  virtual void set_type(const SurfacePropertyType new_type_) {
    if (initialized) {
      throw RuntimeError("Value 'type' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    type = new_type_;
  }
  virtual SurfacePropertyType get_type() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return type;
  }

  std::shared_ptr<Complex> affected_complex_pattern;
  virtual void set_affected_complex_pattern(std::shared_ptr<Complex> new_affected_complex_pattern_) {
    if (initialized) {
      throw RuntimeError("Value 'affected_complex_pattern' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    affected_complex_pattern = new_affected_complex_pattern_;
  }
  virtual std::shared_ptr<Complex> get_affected_complex_pattern() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return affected_complex_pattern;
  }

  float_t concentration;
  virtual void set_concentration(const float_t new_concentration_) {
    if (initialized) {
      throw RuntimeError("Value 'concentration' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    concentration = new_concentration_;
  }
  virtual float_t get_concentration() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return concentration;
  }

  // --- methods ---
}; // GenSurfaceProperty

class SurfaceProperty;
py::class_<SurfaceProperty> define_pybinding_SurfaceProperty(py::module& m);
} // namespace API
} // namespace MCell

#endif // API_GEN_SURFACE_PROPERTY_H
