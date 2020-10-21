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

#ifndef API_GEN_GEOMETRY_OBJECT_H
#define API_GEN_GEOMETRY_OBJECT_H

#include "../api/common.h"
#include "../api/region.h"


namespace MCell {
namespace API {

class InitialSurfaceRelease;
class Region;
class SurfaceClass;
class SurfaceRegion;

#define GEOMETRY_OBJECT_CTOR() \
    GeometryObject( \
        const std::string& name_, \
        const std::vector<std::vector<float_t>> vertex_list_, \
        const std::vector<std::vector<int>> element_connections_, \
        const bool is_bngl_compartment_ = false, \
        const std::string& surface_compartment_name_ = STR_UNSET, \
        const std::vector<std::shared_ptr<SurfaceRegion>> surface_regions_ = std::vector<std::shared_ptr<SurfaceRegion>>(), \
        std::shared_ptr<SurfaceClass> surface_class_ = nullptr, \
        const std::vector<std::shared_ptr<InitialSurfaceRelease>> initial_surface_releases_ = std::vector<std::shared_ptr<InitialSurfaceRelease>>(), \
        const RegionNodeType node_type_ = RegionNodeType::UNSET, \
        std::shared_ptr<Region> left_node_ = nullptr, \
        std::shared_ptr<Region> right_node_ = nullptr \
    )  : GenGeometryObject(node_type_,left_node_,right_node_) { \
      class_name = "GeometryObject"; \
      name = name_; \
      vertex_list = vertex_list_; \
      element_connections = element_connections_; \
      is_bngl_compartment = is_bngl_compartment_; \
      surface_compartment_name = surface_compartment_name_; \
      surface_regions = surface_regions_; \
      surface_class = surface_class_; \
      initial_surface_releases = initial_surface_releases_; \
      node_type = node_type_; \
      left_node = left_node_; \
      right_node = right_node_; \
      postprocess_in_ctor();\
      check_semantics();\
    }

class GenGeometryObject: public Region {
public:
  GenGeometryObject( 
      const RegionNodeType node_type_ = RegionNodeType::UNSET, 
      std::shared_ptr<Region> left_node_ = nullptr, 
      std::shared_ptr<Region> right_node_ = nullptr 
  )  : Region(node_type_,left_node_,right_node_)  {
  }
  void postprocess_in_ctor() override {}
  void check_semantics() const override;
  void set_initialized() override;
  bool __eq__(const GenGeometryObject& other) const;
  void set_all_attributes_as_default_or_unset() override;

  std::string to_str(const std::string ind="") const override;

  // --- attributes ---
  std::vector<std::vector<float_t>> vertex_list;
  virtual void set_vertex_list(const std::vector<std::vector<float_t>> new_vertex_list_) {
    if (initialized) {
      throw RuntimeError("Value 'vertex_list' of object with name " + name + " (class " + class_name + ")"
                         "cannot be set after model was initialized.");
    }
    vertex_list = new_vertex_list_;
  }
  virtual std::vector<std::vector<float_t>> get_vertex_list() const {
    return vertex_list;
  }

  std::vector<std::vector<int>> element_connections;
  virtual void set_element_connections(const std::vector<std::vector<int>> new_element_connections_) {
    if (initialized) {
      throw RuntimeError("Value 'element_connections' of object with name " + name + " (class " + class_name + ")"
                         "cannot be set after model was initialized.");
    }
    element_connections = new_element_connections_;
  }
  virtual std::vector<std::vector<int>> get_element_connections() const {
    return element_connections;
  }

  bool is_bngl_compartment;
  virtual void set_is_bngl_compartment(const bool new_is_bngl_compartment_) {
    if (initialized) {
      throw RuntimeError("Value 'is_bngl_compartment' of object with name " + name + " (class " + class_name + ")"
                         "cannot be set after model was initialized.");
    }
    is_bngl_compartment = new_is_bngl_compartment_;
  }
  virtual bool get_is_bngl_compartment() const {
    return is_bngl_compartment;
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

  std::vector<std::shared_ptr<SurfaceRegion>> surface_regions;
  virtual void set_surface_regions(const std::vector<std::shared_ptr<SurfaceRegion>> new_surface_regions_) {
    if (initialized) {
      throw RuntimeError("Value 'surface_regions' of object with name " + name + " (class " + class_name + ")"
                         "cannot be set after model was initialized.");
    }
    surface_regions = new_surface_regions_;
  }
  virtual std::vector<std::shared_ptr<SurfaceRegion>> get_surface_regions() const {
    return surface_regions;
  }

  std::shared_ptr<SurfaceClass> surface_class;
  virtual void set_surface_class(std::shared_ptr<SurfaceClass> new_surface_class_) {
    if (initialized) {
      throw RuntimeError("Value 'surface_class' of object with name " + name + " (class " + class_name + ")"
                         "cannot be set after model was initialized.");
    }
    surface_class = new_surface_class_;
  }
  virtual std::shared_ptr<SurfaceClass> get_surface_class() const {
    return surface_class;
  }

  std::vector<std::shared_ptr<InitialSurfaceRelease>> initial_surface_releases;
  virtual void set_initial_surface_releases(const std::vector<std::shared_ptr<InitialSurfaceRelease>> new_initial_surface_releases_) {
    if (initialized) {
      throw RuntimeError("Value 'initial_surface_releases' of object with name " + name + " (class " + class_name + ")"
                         "cannot be set after model was initialized.");
    }
    initial_surface_releases = new_initial_surface_releases_;
  }
  virtual std::vector<std::shared_ptr<InitialSurfaceRelease>> get_initial_surface_releases() const {
    return initial_surface_releases;
  }

  // --- methods ---
  virtual void translate(const Vec3& move) = 0;
}; // GenGeometryObject

class GeometryObject;
py::class_<GeometryObject> define_pybinding_GeometryObject(py::module& m);
} // namespace API
} // namespace MCell

#endif // API_GEN_GEOMETRY_OBJECT_H
