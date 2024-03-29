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

#ifndef API_GEN_SURFACE_REGION_H
#define API_GEN_SURFACE_REGION_H

#include "api/api_common.h"
#include "api/region.h"


namespace MCell {
namespace API {

class SurfaceRegion;
class Color;
class InitialSurfaceRelease;
class Region;
class SurfaceClass;
class PythonExportContext;

#define SURFACE_REGION_CTOR() \
    SurfaceRegion( \
        const std::string& name_, \
        const std::vector<int> wall_indices_, \
        std::shared_ptr<SurfaceClass> surface_class_ = nullptr, \
        const std::vector<std::shared_ptr<InitialSurfaceRelease>> initial_surface_releases_ = std::vector<std::shared_ptr<InitialSurfaceRelease>>(), \
        std::shared_ptr<Color> initial_color_ = nullptr, \
        const RegionNodeType node_type_ = RegionNodeType::UNSET, \
        std::shared_ptr<Region> left_node_ = nullptr, \
        std::shared_ptr<Region> right_node_ = nullptr \
    )  : GenSurfaceRegion(node_type_,left_node_,right_node_) { \
      class_name = "SurfaceRegion"; \
      name = name_; \
      wall_indices = wall_indices_; \
      surface_class = surface_class_; \
      initial_surface_releases = initial_surface_releases_; \
      initial_color = initial_color_; \
      node_type = node_type_; \
      left_node = left_node_; \
      right_node = right_node_; \
      postprocess_in_ctor(); \
      check_semantics(); \
    } \
    SurfaceRegion(DefaultCtorArgType) : \
      GenSurfaceRegion(DefaultCtorArgType()) { \
      set_all_attributes_as_default_or_unset(); \
      set_all_custom_attributes_to_default(); \
    }

class GenSurfaceRegion: public Region {
public:
  GenSurfaceRegion( 
      const RegionNodeType node_type_ = RegionNodeType::UNSET, 
      std::shared_ptr<Region> left_node_ = nullptr, 
      std::shared_ptr<Region> right_node_ = nullptr 
  )  : Region(node_type_,left_node_,right_node_)  {
  }
  GenSurfaceRegion(DefaultCtorArgType) {
  }
  void postprocess_in_ctor() override {}
  void check_semantics() const override;
  void set_initialized() override;
  void set_all_attributes_as_default_or_unset() override;

  std::shared_ptr<SurfaceRegion> copy_surface_region() const;
  std::shared_ptr<SurfaceRegion> deepcopy_surface_region(py::dict = py::dict()) const;
  virtual bool __eq__(const SurfaceRegion& other) const;
  virtual bool eq_nonarray_attributes(const SurfaceRegion& other, const bool ignore_name = false) const;
  bool operator == (const SurfaceRegion& other) const { return __eq__(other);}
  bool operator != (const SurfaceRegion& other) const { return !__eq__(other);}
  std::string to_str(const bool all_details=false, const std::string ind="") const override;

  virtual std::string export_to_python(std::ostream& out, PythonExportContext& ctx);
  virtual std::string export_vec_wall_indices(std::ostream& out, PythonExportContext& ctx, const std::string& parent_name);
  virtual std::string export_vec_initial_surface_releases(std::ostream& out, PythonExportContext& ctx, const std::string& parent_name);


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
  virtual std::vector<int>& get_wall_indices() {
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
  virtual std::vector<std::shared_ptr<InitialSurfaceRelease>>& get_initial_surface_releases() {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return initial_surface_releases;
  }

  std::shared_ptr<Color> initial_color;
  virtual void set_initial_color(std::shared_ptr<Color> new_initial_color_) {
    if (initialized) {
      throw RuntimeError("Value 'initial_color' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    initial_color = new_initial_color_;
  }
  virtual std::shared_ptr<Color> get_initial_color() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return initial_color;
  }

  // --- methods ---
}; // GenSurfaceRegion

class SurfaceRegion;
py::class_<SurfaceRegion> define_pybinding_SurfaceRegion(py::module& m);
} // namespace API
} // namespace MCell

#endif // API_GEN_SURFACE_REGION_H
