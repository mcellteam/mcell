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

#ifndef API_GEN_WALL_H
#define API_GEN_WALL_H

#include "api/common.h"
#include "api/base_introspection_class.h"

namespace MCell {
namespace API {

class Wall;
class GeometryObject;

#define WALL_CTOR_NOARGS() \
    Wall( \
    ) { \
      class_name = "Wall"; \
      geometry_object = nullptr; \
      wall_index = INT_UNSET; \
      vertices = std::vector<Vec3>(); \
      area = FLT_UNSET; \
      normal = VEC3_UNSET; \
      is_movable = true; \
      postprocess_in_ctor();\
      check_semantics();\
    }

class GenWall: public BaseIntrospectionClass {
public:
  void postprocess_in_ctor() override {}
  void check_semantics() const override;
  void set_initialized() override;
  void set_all_attributes_as_default_or_unset() override;

  virtual bool __eq__(const Wall& other) const;
  virtual bool eq_nonarray_attributes(const Wall& other) const;
  bool operator == (const Wall& other) const { return __eq__(other);}
  bool operator != (const Wall& other) const { return !__eq__(other);}
  std::string to_str(const std::string ind="") const override;

  // --- attributes ---
  std::shared_ptr<GeometryObject> geometry_object;
  virtual void set_geometry_object(std::shared_ptr<GeometryObject> new_geometry_object_) {
    if (initialized) {
      throw RuntimeError("Value 'geometry_object' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    geometry_object = new_geometry_object_;
  }
  virtual std::shared_ptr<GeometryObject> get_geometry_object() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return geometry_object;
  }

  int wall_index;
  virtual void set_wall_index(const int new_wall_index_) {
    if (initialized) {
      throw RuntimeError("Value 'wall_index' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    wall_index = new_wall_index_;
  }
  virtual int get_wall_index() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return wall_index;
  }

  std::vector<Vec3> vertices;
  virtual void set_vertices(const std::vector<Vec3> new_vertices_) {
    if (initialized) {
      throw RuntimeError("Value 'vertices' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    vertices = new_vertices_;
  }
  virtual std::vector<Vec3> get_vertices() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return vertices;
  }

  float_t area;
  virtual void set_area(const float_t new_area_) {
    if (initialized) {
      throw RuntimeError("Value 'area' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    area = new_area_;
  }
  virtual float_t get_area() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return area;
  }

  Vec3 normal;
  virtual void set_normal(const Vec3& new_normal_) {
    if (initialized) {
      throw RuntimeError("Value 'normal' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    normal = new_normal_;
  }
  virtual const Vec3& get_normal() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return normal;
  }

  bool is_movable;
  virtual void set_is_movable(const bool new_is_movable_) {
    if (initialized) {
      throw RuntimeError("Value 'is_movable' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    is_movable = new_is_movable_;
  }
  virtual bool get_is_movable() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return is_movable;
  }

  // --- methods ---
}; // GenWall

class Wall;
py::class_<Wall> define_pybinding_Wall(py::module& m);
} // namespace API
} // namespace MCell

#endif // API_GEN_WALL_H
