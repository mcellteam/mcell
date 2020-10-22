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

#include "../api/common.h"
#include "../api/base_data_class.h"

namespace MCell {
namespace API {

class Species;

#define SURFACE_PROPERTY_CTOR() \
    SurfaceProperty( \
        const SurfacePropertyType type_ = SurfacePropertyType::UNSET, \
        std::shared_ptr<Species> affected_species_ = nullptr, \
        const Orientation orientation_ = Orientation::NOT_SET \
    ) { \
      class_name = "SurfaceProperty"; \
      type = type_; \
      affected_species = affected_species_; \
      orientation = orientation_; \
      postprocess_in_ctor();\
      check_semantics();\
    }

class GenSurfaceProperty: public BaseDataClass {
public:
  void postprocess_in_ctor() override {}
  void check_semantics() const override;
  void set_initialized() override;
  bool __eq__(const GenSurfaceProperty& other) const;
  void set_all_attributes_as_default_or_unset() override;

  std::string to_str(const std::string ind="") const override;

  // --- attributes ---
  SurfacePropertyType type;
  virtual void set_type(const SurfacePropertyType new_type_) {
    if (initialized) {
      throw RuntimeError("Value 'type' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    type = new_type_;
  }
  virtual SurfacePropertyType get_type() const {
    return type;
  }

  std::shared_ptr<Species> affected_species;
  virtual void set_affected_species(std::shared_ptr<Species> new_affected_species_) {
    if (initialized) {
      throw RuntimeError("Value 'affected_species' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    affected_species = new_affected_species_;
  }
  virtual std::shared_ptr<Species> get_affected_species() const {
    return affected_species;
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

  // --- methods ---
}; // GenSurfaceProperty

class SurfaceProperty;
py::class_<SurfaceProperty> define_pybinding_SurfaceProperty(py::module& m);
} // namespace API
} // namespace MCell

#endif // API_GEN_SURFACE_PROPERTY_H
