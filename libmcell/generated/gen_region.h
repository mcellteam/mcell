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

#ifndef API_GEN_REGION_H
#define API_GEN_REGION_H

#include "../api/common.h"
#include "../api/base_data_class.h"

namespace MCell {
namespace API {

class GeometryObject;
class Region;
class SurfaceArea;

#define REGION_CTOR() \
    Region( \
        const RegionNodeType node_type_ = RegionNodeType::Unset, \
        std::shared_ptr<GeometryObject> geometry_object_ = nullptr, \
        std::shared_ptr<SurfaceArea> surface_area_ = nullptr, \
        std::shared_ptr<Region> left_node_ = nullptr, \
        std::shared_ptr<Region> right_node_ = nullptr \
    ) { \
      class_name = "Region"; \
      node_type = node_type_; \
      geometry_object = geometry_object_; \
      surface_area = surface_area_; \
      left_node = left_node_; \
      right_node = right_node_; \
      postprocess_in_ctor();\
      check_semantics();\
    }

class GenRegion: public BaseDataClass {
public:
  void postprocess_in_ctor() override {}
  void check_semantics() const override;
  void set_initialized() override;

  bool __eq__(const GenRegion& other) const;
  std::string to_str(const std::string ind="") const override;

  // --- attributes ---
  RegionNodeType node_type;
  virtual void set_node_type(const RegionNodeType new_node_type_) {
    if (initialized) {
      throw RuntimeError("Value 'node_type' of object with name " + name + " (class " + class_name + ")"
                         "cannot be set after model was initialized.");
    }
    node_type = new_node_type_;
  }
  virtual RegionNodeType get_node_type() const {
    return node_type;
  }

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

  std::shared_ptr<SurfaceArea> surface_area;
  virtual void set_surface_area(std::shared_ptr<SurfaceArea> new_surface_area_) {
    if (initialized) {
      throw RuntimeError("Value 'surface_area' of object with name " + name + " (class " + class_name + ")"
                         "cannot be set after model was initialized.");
    }
    surface_area = new_surface_area_;
  }
  virtual std::shared_ptr<SurfaceArea> get_surface_area() const {
    return surface_area;
  }

  std::shared_ptr<Region> left_node;
  virtual void set_left_node(std::shared_ptr<Region> new_left_node_) {
    if (initialized) {
      throw RuntimeError("Value 'left_node' of object with name " + name + " (class " + class_name + ")"
                         "cannot be set after model was initialized.");
    }
    left_node = new_left_node_;
  }
  virtual std::shared_ptr<Region> get_left_node() const {
    return left_node;
  }

  std::shared_ptr<Region> right_node;
  virtual void set_right_node(std::shared_ptr<Region> new_right_node_) {
    if (initialized) {
      throw RuntimeError("Value 'right_node' of object with name " + name + " (class " + class_name + ")"
                         "cannot be set after model was initialized.");
    }
    right_node = new_right_node_;
  }
  virtual std::shared_ptr<Region> get_right_node() const {
    return right_node;
  }

  // --- methods ---
  virtual std::shared_ptr<Region> union_(std::shared_ptr<Region> second_region) = 0;
  virtual std::shared_ptr<Region> difference(std::shared_ptr<Region> second_region) = 0;
  virtual std::shared_ptr<Region> intersect(std::shared_ptr<Region> second_region) = 0;
}; // GenRegion

class Region;
py::class_<Region> define_pybinding_Region(py::module& m);
} // namespace API
} // namespace MCell

#endif // API_GEN_REGION_H
