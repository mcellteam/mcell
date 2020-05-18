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

namespace MCell {
namespace API {

class SurfaceRegion;

#define GEOMETRY_OBJECT_CTOR() \
    GeometryObject( \
        const std::string& name_, \
        const std::vector<std::vector<float_t>> vertex_list_, \
        const std::vector<std::vector<int>> element_connections_, \
        const std::vector<std::shared_ptr<SurfaceRegion>> surface_regions_ = std::vector<std::shared_ptr<SurfaceRegion>>() \
    ) { \
      class_name = "GeometryObject"; \
      name = name_; \
      vertex_list = vertex_list_; \
      element_connections = element_connections_; \
      surface_regions = surface_regions_; \
      postprocess_in_ctor();\
      check_semantics();\
    }

class GenGeometryObject: public BaseDataClass {
public:
  void postprocess_in_ctor() override {}
  void check_semantics() const override;
  std::string to_str(const std::string ind="") const override;

  // --- attributes ---
  std::vector<std::vector<float_t>> vertex_list;
  virtual void set_vertex_list(const std::vector<std::vector<float_t>> new_vertex_list_) {
    vertex_list = new_vertex_list_;
  }
  virtual std::vector<std::vector<float_t>> get_vertex_list() const {
    return vertex_list;
  }

  std::vector<std::vector<int>> element_connections;
  virtual void set_element_connections(const std::vector<std::vector<int>> new_element_connections_) {
    element_connections = new_element_connections_;
  }
  virtual std::vector<std::vector<int>> get_element_connections() const {
    return element_connections;
  }

  std::vector<std::shared_ptr<SurfaceRegion>> surface_regions;
  virtual void set_surface_regions(const std::vector<std::shared_ptr<SurfaceRegion>> new_surface_regions_) {
    surface_regions = new_surface_regions_;
  }
  virtual std::vector<std::shared_ptr<SurfaceRegion>> get_surface_regions() const {
    return surface_regions;
  }

  // --- methods ---
}; // GenGeometryObject

class GeometryObject;
py::class_<GeometryObject> define_pybinding_GeometryObject(py::module& m);
} // namespace API
} // namespace MCell

#endif // API_GEN_GEOMETRY_OBJECT_H
