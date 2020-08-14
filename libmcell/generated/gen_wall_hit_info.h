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

#ifndef API_GEN_WALL_HIT_INFO_H
#define API_GEN_WALL_HIT_INFO_H

#include "../api/common.h"

namespace MCell {
namespace API {

class GenWallHitInfo {
public:
  virtual ~GenWallHitInfo() {}
  std::string to_str(const std::string ind="") const ;

  // --- attributes ---
  int molecule_id;
  virtual void set_molecule_id(const int new_molecule_id_) {
    molecule_id = new_molecule_id_;
  }
  virtual int get_molecule_id() const {
    return molecule_id;
  }

  int geometry_object_id;
  virtual void set_geometry_object_id(const int new_geometry_object_id_) {
    geometry_object_id = new_geometry_object_id_;
  }
  virtual int get_geometry_object_id() const {
    return geometry_object_id;
  }

  int wall_id;
  virtual void set_wall_id(const int new_wall_id_) {
    wall_id = new_wall_id_;
  }
  virtual int get_wall_id() const {
    return wall_id;
  }

  float_t time;
  virtual void set_time(const float_t new_time_) {
    time = new_time_;
  }
  virtual float_t get_time() const {
    return time;
  }

  Vec3 pos;
  virtual void set_pos(const Vec3& new_pos_) {
    pos = new_pos_;
  }
  virtual const Vec3& get_pos() const {
    return pos;
  }

  float_t time_before_hit;
  virtual void set_time_before_hit(const float_t new_time_before_hit_) {
    time_before_hit = new_time_before_hit_;
  }
  virtual float_t get_time_before_hit() const {
    return time_before_hit;
  }

  Vec3 pos_before_hit;
  virtual void set_pos_before_hit(const Vec3& new_pos_before_hit_) {
    pos_before_hit = new_pos_before_hit_;
  }
  virtual const Vec3& get_pos_before_hit() const {
    return pos_before_hit;
  }

  // --- methods ---
}; // GenWallHitInfo

class WallHitInfo;
py::class_<WallHitInfo> define_pybinding_WallHitInfo(py::module& m);
} // namespace API
} // namespace MCell

#endif // API_GEN_WALL_HIT_INFO_H
