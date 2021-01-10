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

#ifndef API_GEN_INITIAL_SURFACE_RELEASE_H
#define API_GEN_INITIAL_SURFACE_RELEASE_H

#include "api/common.h"
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
        const float_t density_ = FLT_UNSET \
    ) { \
      class_name = "InitialSurfaceRelease"; \
      complex = complex_; \
      number_to_release = number_to_release_; \
      density = density_; \
      postprocess_in_ctor();\
      check_semantics();\
    }

class GenInitialSurfaceRelease: public BaseDataClass {
public:
  void postprocess_in_ctor() override {}
  void check_semantics() const override;
  void set_initialized() override;
  void set_all_attributes_as_default_or_unset() override;

  virtual bool __eq__(const InitialSurfaceRelease& other) const;
  virtual bool eq_nonarray_attributes(const InitialSurfaceRelease& other, const bool ignore_name = false) const;
  bool operator == (const InitialSurfaceRelease& other) const { return __eq__(other);}
  bool operator != (const InitialSurfaceRelease& other) const { return !__eq__(other);}
  std::string to_str(const std::string ind="") const override;

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

  float_t density;
  virtual void set_density(const float_t new_density_) {
    if (initialized) {
      throw RuntimeError("Value 'density' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    density = new_density_;
  }
  virtual float_t get_density() const {
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
