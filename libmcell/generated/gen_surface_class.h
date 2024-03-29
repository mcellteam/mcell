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

#ifndef API_GEN_SURFACE_CLASS_H
#define API_GEN_SURFACE_CLASS_H

#include "api/api_common.h"
#include "api/surface_property.h"


namespace MCell {
namespace API {

class SurfaceClass;
class Complex;
class SurfaceProperty;
class PythonExportContext;

#define SURFACE_CLASS_CTOR() \
    SurfaceClass( \
        const std::string& name_, \
        const std::vector<std::shared_ptr<SurfaceProperty>> properties_ = std::vector<std::shared_ptr<SurfaceProperty>>(), \
        const SurfacePropertyType type_ = SurfacePropertyType::UNSET, \
        std::shared_ptr<Complex> affected_complex_pattern_ = nullptr, \
        const double concentration_ = FLT_UNSET \
    )  : GenSurfaceClass(type_,affected_complex_pattern_,concentration_) { \
      class_name = "SurfaceClass"; \
      name = name_; \
      properties = properties_; \
      type = type_; \
      affected_complex_pattern = affected_complex_pattern_; \
      concentration = concentration_; \
      postprocess_in_ctor(); \
      check_semantics(); \
    } \
    SurfaceClass(DefaultCtorArgType) : \
      GenSurfaceClass(DefaultCtorArgType()) { \
      set_all_attributes_as_default_or_unset(); \
      set_all_custom_attributes_to_default(); \
    }

class GenSurfaceClass: public SurfaceProperty {
public:
  GenSurfaceClass( 
      const SurfacePropertyType type_ = SurfacePropertyType::UNSET, 
      std::shared_ptr<Complex> affected_complex_pattern_ = nullptr, 
      const double concentration_ = FLT_UNSET 
  )  : SurfaceProperty(type_,affected_complex_pattern_,concentration_)  {
  }
  GenSurfaceClass(DefaultCtorArgType) {
  }
  void postprocess_in_ctor() override {}
  void check_semantics() const override;
  void set_initialized() override;
  void set_all_attributes_as_default_or_unset() override;

  std::shared_ptr<SurfaceClass> copy_surface_class() const;
  std::shared_ptr<SurfaceClass> deepcopy_surface_class(py::dict = py::dict()) const;
  virtual bool __eq__(const SurfaceClass& other) const;
  virtual bool eq_nonarray_attributes(const SurfaceClass& other, const bool ignore_name = false) const;
  bool operator == (const SurfaceClass& other) const { return __eq__(other);}
  bool operator != (const SurfaceClass& other) const { return !__eq__(other);}
  std::string to_str(const bool all_details=false, const std::string ind="") const override;

  virtual std::string export_to_python(std::ostream& out, PythonExportContext& ctx);
  virtual std::string export_vec_properties(std::ostream& out, PythonExportContext& ctx, const std::string& parent_name);


  // --- attributes ---
  std::vector<std::shared_ptr<SurfaceProperty>> properties;
  virtual void set_properties(const std::vector<std::shared_ptr<SurfaceProperty>> new_properties_) {
    if (initialized) {
      throw RuntimeError("Value 'properties' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    properties = new_properties_;
  }
  virtual std::vector<std::shared_ptr<SurfaceProperty>>& get_properties() {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return properties;
  }

  // --- methods ---
}; // GenSurfaceClass

class SurfaceClass;
py::class_<SurfaceClass> define_pybinding_SurfaceClass(py::module& m);
} // namespace API
} // namespace MCell

#endif // API_GEN_SURFACE_CLASS_H
