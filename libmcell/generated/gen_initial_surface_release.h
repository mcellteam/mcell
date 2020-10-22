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

#include "../api/common.h"
#include "../api/base_data_class.h"

namespace MCell {
namespace API {

class Species;

#define INITIAL_SURFACE_RELEASE_CTOR() \
    InitialSurfaceRelease( \
        std::shared_ptr<Species> species_ = nullptr, \
        const std::string& bngl_species_ = STR_UNSET, \
        const Orientation orientation_ = Orientation::UP, \
        const int number_to_release_ = INT_UNSET, \
        const float_t density_ = FLT_UNSET \
    ) { \
      class_name = "InitialSurfaceRelease"; \
      species = species_; \
      bngl_species = bngl_species_; \
      orientation = orientation_; \
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
  bool __eq__(const GenInitialSurfaceRelease& other) const;
  void set_all_attributes_as_default_or_unset() override;

  std::string to_str(const std::string ind="") const override;

  // --- attributes ---
  std::shared_ptr<Species> species;
  virtual void set_species(std::shared_ptr<Species> new_species_) {
    if (initialized) {
      throw RuntimeError("Value 'species' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    species = new_species_;
  }
  virtual std::shared_ptr<Species> get_species() const {
    return species;
  }

  std::string bngl_species;
  virtual void set_bngl_species(const std::string& new_bngl_species_) {
    if (initialized) {
      throw RuntimeError("Value 'bngl_species' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    bngl_species = new_bngl_species_;
  }
  virtual const std::string& get_bngl_species() const {
    return bngl_species;
  }

  Orientation orientation;
  virtual void set_orientation(const Orientation new_orientation_) {
    if (initialized) {
      throw RuntimeError("Value 'orientation' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    orientation = new_orientation_;
  }
  virtual Orientation get_orientation() const {
    return orientation;
  }

  int number_to_release;
  virtual void set_number_to_release(const int new_number_to_release_) {
    if (initialized) {
      throw RuntimeError("Value 'number_to_release' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    number_to_release = new_number_to_release_;
  }
  virtual int get_number_to_release() const {
    return number_to_release;
  }

  float_t density;
  virtual void set_density(const float_t new_density_) {
    if (initialized) {
      throw RuntimeError("Value 'density' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    density = new_density_;
  }
  virtual float_t get_density() const {
    return density;
  }

  // --- methods ---
}; // GenInitialSurfaceRelease

class InitialSurfaceRelease;
py::class_<InitialSurfaceRelease> define_pybinding_InitialSurfaceRelease(py::module& m);
} // namespace API
} // namespace MCell

#endif // API_GEN_INITIAL_SURFACE_RELEASE_H
