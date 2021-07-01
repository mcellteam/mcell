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

#ifndef API_GEN_INITIAL_SURFACE_RELEASE_H
#define API_GEN_INITIAL_SURFACE_RELEASE_H

#include "api/api_common.h"
#include "api/base_data_class.h"

namespace MCell {
namespace API {

class InitialSurfaceRelease;
class Complex;
class PythonExportContext;

#define INITIAL_SURFACE_RELEASE_CTOR() \
    InitialSurfaceRelease( \
        std::shared_ptr<Complex> complex_, \
        const int number_to_release_ = INT_UNSET, \
        const double density_ = FLT_UNSET \
    ) { \
      class_name = "InitialSurfaceRelease"; \
      complex = complex_; \
      number_to_release = number_to_release_; \
      density = density_; \
      postprocess_in_ctor(); \
      check_semantics(); \
    } \
    InitialSurfaceRelease(DefaultCtorArgType) : \
      GenInitialSurfaceRelease(DefaultCtorArgType()) { \
      set_all_attributes_as_default_or_unset(); \
      set_all_custom_attributes_to_default(); \
    }

class GenInitialSurfaceRelease: public BaseDataClass {
public:
  GenInitialSurfaceRelease() {
  }
  GenInitialSurfaceRelease(DefaultCtorArgType) {
  }
  void postprocess_in_ctor() override {}
  void check_semantics() const override;
  void set_initialized() override;
  void set_all_attributes_as_default_or_unset() override;

  std::shared_ptr<InitialSurfaceRelease> copy_initial_surface_release() const;
  std::shared_ptr<InitialSurfaceRelease> deepcopy_initial_surface_release(py::dict = py::dict()) const;
  virtual bool __eq__(const InitialSurfaceRelease& other) const;
  virtual bool eq_nonarray_attributes(const InitialSurfaceRelease& other, const bool ignore_name = false) const;
  bool operator == (const InitialSurfaceRelease& other) const { return __eq__(other);}
  bool operator != (const InitialSurfaceRelease& other) const { return !__eq__(other);}
  std::string to_str(const bool all_details=false, const std::string ind="") const override;

  std::string export_to_python(std::ostream& out, PythonExportContext& ctx) override;


  // --- attributes ---
  std::shared_ptr<Complex> complex;
  virtual void set_complex(std::shared_ptr<Complex> new_complex_) {
    if (initialized) {
      throw RuntimeError("Value 'complex' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    complex = new_complex_;
  }
  virtual std::shared_ptr<Complex> get_complex() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return complex;
  }

  int number_to_release;
  virtual void set_number_to_release(const int new_number_to_release_) {
    if (initialized) {
      throw RuntimeError("Value 'number_to_release' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    number_to_release = new_number_to_release_;
  }
  virtual int get_number_to_release() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return number_to_release;
  }

  double density;
  virtual void set_density(const double new_density_) {
    if (initialized) {
      throw RuntimeError("Value 'density' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    density = new_density_;
  }
  virtual double get_density() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return density;
  }

  // --- methods ---
}; // GenInitialSurfaceRelease

class InitialSurfaceRelease;
py::class_<InitialSurfaceRelease> define_pybinding_InitialSurfaceRelease(py::module& m);
} // namespace API
} // namespace MCell

#endif // API_GEN_INITIAL_SURFACE_RELEASE_H
