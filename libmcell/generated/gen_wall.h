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

#ifndef API_GEN_WALL_H
#define API_GEN_WALL_H

#include "api/api_common.h"
#include "api/base_introspection_class.h"

namespace MCell {
namespace API {

class Wall;
class GeometryObject;
class PythonExportContext;

#define WALL_CTOR_NOARGS() \
    Wall( \
    ) { \
      class_name = "Wall"; \
      geometry_object = nullptr; \
      wall_index = INT_UNSET; \
      vertices = std::vector<std::vector<double>>(); \
      area = FLT_UNSET; \
      unit_normal = std::vector<double>(); \
      is_movable = true; \
      postprocess_in_ctor(); \
      check_semantics(); \
    } \
    Wall(DefaultCtorArgType) : \
      GenWall(DefaultCtorArgType()) { \
      set_all_attributes_as_default_or_unset(); \
      set_all_custom_attributes_to_default(); \
    }

class GenWall: public BaseIntrospectionClass {
public:
  GenWall() {
  }
  GenWall(DefaultCtorArgType) {
  }
  void postprocess_in_ctor() override {}
  void check_semantics() const override;
  void set_initialized() override;
  void set_all_attributes_as_default_or_unset() override;

  std::shared_ptr<Wall> copy_wall() const;
  std::shared_ptr<Wall> deepcopy_wall(py::dict = py::dict()) const;
  virtual bool __eq__(const Wall& other) const;
  virtual bool eq_nonarray_attributes(const Wall& other, const bool ignore_name = false) const;
  bool operator == (const Wall& other) const { return __eq__(other);}
  bool operator != (const Wall& other) const { return !__eq__(other);}
  std::string to_str(const bool all_details=false, const std::string ind="") const override;

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

  std::vector<std::vector<double>> vertices;
  virtual void set_vertices(const std::vector<std::vector<double>> new_vertices_) {
    if (initialized) {
      throw RuntimeError("Value 'vertices' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    vertices = new_vertices_;
  }
  virtual std::vector<std::vector<double>>& get_vertices() {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return vertices;
  }

  double area;
  virtual void set_area(const double new_area_) {
    if (initialized) {
      throw RuntimeError("Value 'area' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    area = new_area_;
  }
  virtual double get_area() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return area;
  }

  std::vector<double> unit_normal;
  virtual void set_unit_normal(const std::vector<double> new_unit_normal_) {
    if (initialized) {
      throw RuntimeError("Value 'unit_normal' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    unit_normal = new_unit_normal_;
  }
  virtual std::vector<double>& get_unit_normal() {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return unit_normal;
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
