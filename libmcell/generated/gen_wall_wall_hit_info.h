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

#ifndef API_GEN_WALL_WALL_HIT_INFO_H
#define API_GEN_WALL_WALL_HIT_INFO_H

#include "../api/common.h"
#include "../api/base_introspection_class.h"

namespace MCell {
namespace API {

class Wall;

#define WALL_WALL_HIT_INFO_CTOR_NOARGS() \
    WallWallHitInfo( \
    ) { \
      class_name = "WallWallHitInfo"; \
      wall1 = nullptr; \
      wall2 = nullptr; \
      postprocess_in_ctor();\
      check_semantics();\
    }

class GenWallWallHitInfo: public BaseIntrospectionClass {
public:
  void postprocess_in_ctor() override {}
  void check_semantics() const override;
  void set_initialized() override;
  bool __eq__(const GenWallWallHitInfo& other) const;
  void set_all_attributes_as_default_or_unset() override;

  std::string to_str(const std::string ind="") const override;

  // --- attributes ---
  std::shared_ptr<Wall> wall1;
  virtual void set_wall1(std::shared_ptr<Wall> new_wall1_) {
    if (initialized) {
      throw RuntimeError("Value 'wall1' of object with name " + name + " (class " + class_name + ")"
                         "cannot be set after model was initialized.");
    }
    wall1 = new_wall1_;
  }
  virtual std::shared_ptr<Wall> get_wall1() const {
    return wall1;
  }

  std::shared_ptr<Wall> wall2;
  virtual void set_wall2(std::shared_ptr<Wall> new_wall2_) {
    if (initialized) {
      throw RuntimeError("Value 'wall2' of object with name " + name + " (class " + class_name + ")"
                         "cannot be set after model was initialized.");
    }
    wall2 = new_wall2_;
  }
  virtual std::shared_ptr<Wall> get_wall2() const {
    return wall2;
  }

  // --- methods ---
}; // GenWallWallHitInfo

class WallWallHitInfo;
py::class_<WallWallHitInfo> define_pybinding_WallWallHitInfo(py::module& m);
} // namespace API
} // namespace MCell

#endif // API_GEN_WALL_WALL_HIT_INFO_H
