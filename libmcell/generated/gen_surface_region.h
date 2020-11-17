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

#ifndef API_GEN_SURFACE_REGION_H
#define API_GEN_SURFACE_REGION_H

#include "api/common.h"
#include "api/region.h"


namespace MCell {
namespace API {

class SurfaceRegion;
class InitialSurfaceRelease;
class Region;
class SurfaceClass;

#define SURFACE_REGION_CTOR() \
    SurfaceRegion( \
        const std::string& name_, \
        const std::vector<int> wall_indices_, \
        std::shared_ptr<SurfaceClass> surface_class_ = nullptr, \
        const std::vector<std::shared_ptr<InitialSurfaceRelease>> initial_surface_releases_ = std::vector<std::shared_ptr<InitialSurfaceRelease>>(), \
        const RegionNodeType node_type_ = RegionNodeType::UNSET, \
        std::shared_ptr<Region> left_node_ = nullptr, \
        std::shared_ptr<Region> right_node_ = nullptr \
    )  : GenSurfaceRegion(node_type_,left_node_,right_node_) { \
      class_name = "SurfaceRegion"; \
      name = name_; \
      wall_indices = wall_indices_; \
      surface_class = surface_class_; \
      initial_surface_releases = initial_surface_releases_; \
      node_type = node_type_; \
      left_node = left_node_; \
      right_node = right_node_; \
      postprocess_in_ctor();\
      check_semantics();\
    }

class GenSurfaceRegion: public Region {
public:
  GenSurfaceRegion( 
      const RegionNodeType node_type_ = RegionNodeType::UNSET, 
      std::shared_ptr<Region> left_node_ = nullptr, 
      std::shared_ptr<Region> right_node_ = nullptr 
  )  : Region(node_type_,left_node_,right_node_)  {
  }
  void postprocess_in_ctor() override {}
  void check_semantics() const override;
  void set_initialized() override;
  void set_all_attributes_as_default_or_unset() override;

  virtual bool __eq__(const SurfaceRegion& other) const;
  bool operator == (const SurfaceRegion& other) const { return __eq__(other);}
  bool operator != (const SurfaceRegion& other) const { return !__eq__(other);}
  std::string to_str(const std::string ind="") const override;

  // --- attributes ---
  std::vector<int> wall_indices;
  virtual void set_wall_indices(const std::vector<int> new_wall_indices_) {
    if (initialized) {
      throw RuntimeError("Value 'wall_indices' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    wall_indices = new_wall_indices_;
  }
  virtual std::vector<int> get_wall_indices() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return wall_indices;
  }

  std::shared_ptr<SurfaceClass> surface_class;
  virtual void set_surface_class(std::shared_ptr<SurfaceClass> new_surface_class_) {
    if (initialized) {
      throw RuntimeError("Value 'surface_class' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    surface_class = new_surface_class_;
  }
  virtual std::shared_ptr<SurfaceClass> get_surface_class() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return surface_class;
  }

  std::vector<std::shared_ptr<InitialSurfaceRelease>> initial_surface_releases;
  virtual void set_initial_surface_releases(const std::vector<std::shared_ptr<InitialSurfaceRelease>> new_initial_surface_releases_) {
    if (initialized) {
      throw RuntimeError("Value 'initial_surface_releases' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    initial_surface_releases = new_initial_surface_releases_;
  }
  virtual std::vector<std::shared_ptr<InitialSurfaceRelease>> get_initial_surface_releases() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return initial_surface_releases;
  }

  // --- methods ---
}; // GenSurfaceRegion

class SurfaceRegion;
py::class_<SurfaceRegion> define_pybinding_SurfaceRegion(py::module& m);
} // namespace API
} // namespace MCell

#endif // API_GEN_SURFACE_REGION_H
