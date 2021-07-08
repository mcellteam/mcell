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

#ifndef API_GEN_WALL_WALL_HIT_INFO_H
#define API_GEN_WALL_WALL_HIT_INFO_H

#include "api/api_common.h"
#include "api/base_introspection_class.h"

namespace MCell {
namespace API {

class WallWallHitInfo;
class Wall;
class PythonExportContext;

#define WALL_WALL_HIT_INFO_CTOR_NOARGS() \
    WallWallHitInfo( \
    ) { \
      class_name = "WallWallHitInfo"; \
      wall1 = nullptr; \
      wall2 = nullptr; \
      postprocess_in_ctor(); \
      check_semantics(); \
    } \
    WallWallHitInfo(DefaultCtorArgType) : \
      GenWallWallHitInfo(DefaultCtorArgType()) { \
      set_all_attributes_as_default_or_unset(); \
      set_all_custom_attributes_to_default(); \
    }

class GenWallWallHitInfo: public BaseIntrospectionClass {
public:
  GenWallWallHitInfo() {
  }
  GenWallWallHitInfo(DefaultCtorArgType) {
  }
  void postprocess_in_ctor() override {}
  void check_semantics() const override;
  void set_initialized() override;
  void set_all_attributes_as_default_or_unset() override;

  std::shared_ptr<WallWallHitInfo> copy_wall_wall_hit_info() const;
  std::shared_ptr<WallWallHitInfo> deepcopy_wall_wall_hit_info(py::dict = py::dict()) const;
  virtual bool __eq__(const WallWallHitInfo& other) const;
  virtual bool eq_nonarray_attributes(const WallWallHitInfo& other, const bool ignore_name = false) const;
  bool operator == (const WallWallHitInfo& other) const { return __eq__(other);}
  bool operator != (const WallWallHitInfo& other) const { return !__eq__(other);}
  std::string to_str(const bool all_details=false, const std::string ind="") const override;

  // --- attributes ---
  std::shared_ptr<Wall> wall1;
  virtual void set_wall1(std::shared_ptr<Wall> new_wall1_) {
    if (initialized) {
      throw RuntimeError("Value 'wall1' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    wall1 = new_wall1_;
  }
  virtual std::shared_ptr<Wall> get_wall1() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return wall1;
  }

  std::shared_ptr<Wall> wall2;
  virtual void set_wall2(std::shared_ptr<Wall> new_wall2_) {
    if (initialized) {
      throw RuntimeError("Value 'wall2' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    wall2 = new_wall2_;
  }
  virtual std::shared_ptr<Wall> get_wall2() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return wall2;
  }

  // --- methods ---
}; // GenWallWallHitInfo

class WallWallHitInfo;
py::class_<WallWallHitInfo> define_pybinding_WallWallHitInfo(py::module& m);
} // namespace API
} // namespace MCell

#endif // API_GEN_WALL_WALL_HIT_INFO_H
