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

#ifndef API_GEN_VOLUME_COMPARTMENT_H
#define API_GEN_VOLUME_COMPARTMENT_H

#include "../api/common.h"
#include "../api/base_data_class.h"

namespace MCell {
namespace API {

class GeometryObject;
class Region;
class VolumeCompartment;

#define VOLUME_COMPARTMENT_CTOR() \
    VolumeCompartment( \
        const std::string& name_, \
        std::shared_ptr<GeometryObject> geometry_object_, \
        const std::vector<std::shared_ptr<VolumeCompartment>> child_compartments_ = std::vector<std::shared_ptr<VolumeCompartment>>(), \
        const std::string& surface_compartment_name_ = STR_UNSET \
    ) { \
      class_name = "VolumeCompartment"; \
      name = name_; \
      geometry_object = geometry_object_; \
      child_compartments = child_compartments_; \
      surface_compartment_name = surface_compartment_name_; \
      postprocess_in_ctor();\
      check_semantics();\
    }

class GenVolumeCompartment: public BaseDataClass {
public:
  void postprocess_in_ctor() override {}
  void check_semantics() const override;
  void set_initialized() override;
  bool __eq__(const GenVolumeCompartment& other) const;
  void set_all_attributes_as_default_or_unset() override;

  std::string to_str(const std::string ind="") const override;

  // --- attributes ---
  std::shared_ptr<GeometryObject> geometry_object;
  virtual void set_geometry_object(std::shared_ptr<GeometryObject> new_geometry_object_) {
    if (initialized) {
      throw RuntimeError("Value 'geometry_object' of object with name " + name + " (class " + class_name + ")"
                         "cannot be set after model was initialized.");
    }
    geometry_object = new_geometry_object_;
  }
  virtual std::shared_ptr<GeometryObject> get_geometry_object() const {
    return geometry_object;
  }

  std::vector<std::shared_ptr<VolumeCompartment>> child_compartments;
  virtual void set_child_compartments(const std::vector<std::shared_ptr<VolumeCompartment>> new_child_compartments_) {
    if (initialized) {
      throw RuntimeError("Value 'child_compartments' of object with name " + name + " (class " + class_name + ")"
                         "cannot be set after model was initialized.");
    }
    child_compartments = new_child_compartments_;
  }
  virtual std::vector<std::shared_ptr<VolumeCompartment>> get_child_compartments() const {
    return child_compartments;
  }

  std::string surface_compartment_name;
  virtual void set_surface_compartment_name(const std::string& new_surface_compartment_name_) {
    if (initialized) {
      throw RuntimeError("Value 'surface_compartment_name' of object with name " + name + " (class " + class_name + ")"
                         "cannot be set after model was initialized.");
    }
    surface_compartment_name = new_surface_compartment_name_;
  }
  virtual const std::string& get_surface_compartment_name() const {
    return surface_compartment_name;
  }

  // --- methods ---
  virtual std::shared_ptr<Region> get_volume_compartment_region() = 0;
}; // GenVolumeCompartment

class VolumeCompartment;
py::class_<VolumeCompartment> define_pybinding_VolumeCompartment(py::module& m);
} // namespace API
} // namespace MCell

#endif // API_GEN_VOLUME_COMPARTMENT_H
