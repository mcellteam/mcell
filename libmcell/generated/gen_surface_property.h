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

#ifndef API_GEN_SURFACE_PROPERTY_H
#define API_GEN_SURFACE_PROPERTY_H

#include "api/api_common.h"
#include "api/base_data_class.h"

namespace MCell {
namespace API {

class SurfaceProperty;
class Complex;
class PythonExportContext;

#define SURFACE_PROPERTY_CTOR() \
    SurfaceProperty( \
        const SurfacePropertyType type_ = SurfacePropertyType::UNSET, \
        std::shared_ptr<Complex> affected_complex_pattern_ = nullptr, \
        const double concentration_ = FLT_UNSET \
    ) { \
      class_name = "SurfaceProperty"; \
      type = type_; \
      affected_complex_pattern = affected_complex_pattern_; \
      concentration = concentration_; \
      postprocess_in_ctor(); \
      check_semantics(); \
    } \
    SurfaceProperty(DefaultCtorArgType) : \
      GenSurfaceProperty(DefaultCtorArgType()) { \
      set_all_attributes_as_default_or_unset(); \
      set_all_custom_attributes_to_default(); \
    }

class GenSurfaceProperty: public BaseDataClass {
public:
  GenSurfaceProperty() {
  }
  GenSurfaceProperty(DefaultCtorArgType) {
  }
  void postprocess_in_ctor() override {}
  void check_semantics() const override;
  void set_initialized() override;
  void set_all_attributes_as_default_or_unset() override;

  std::shared_ptr<SurfaceProperty> copy_surface_property() const;
  std::shared_ptr<SurfaceProperty> deepcopy_surface_property(py::dict = py::dict()) const;
  virtual bool __eq__(const SurfaceProperty& other) const;
  virtual bool eq_nonarray_attributes(const SurfaceProperty& other, const bool ignore_name = false) const;
  bool operator == (const SurfaceProperty& other) const { return __eq__(other);}
  bool operator != (const SurfaceProperty& other) const { return !__eq__(other);}
  std::string to_str(const bool all_details=false, const std::string ind="") const override;

  std::string export_to_python(std::ostream& out, PythonExportContext& ctx) override;


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

  double concentration;
  virtual void set_concentration(const double new_concentration_) {
    if (initialized) {
      throw RuntimeError("Value 'concentration' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    concentration = new_concentration_;
  }
  virtual double get_concentration() const {
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
