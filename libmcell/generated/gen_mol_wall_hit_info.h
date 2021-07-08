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

#ifndef API_GEN_MOL_WALL_HIT_INFO_H
#define API_GEN_MOL_WALL_HIT_INFO_H

#include "api/api_common.h"

namespace MCell {
namespace API {

class MolWallHitInfo;
class GeometryObject;
class PythonExportContext;

class GenMolWallHitInfo {
public:
  GenMolWallHitInfo() {
  }
  GenMolWallHitInfo(DefaultCtorArgType) {
  }
  virtual ~GenMolWallHitInfo() {}
  std::shared_ptr<MolWallHitInfo> copy_mol_wall_hit_info() const;
  std::shared_ptr<MolWallHitInfo> deepcopy_mol_wall_hit_info(py::dict = py::dict()) const;
  virtual bool __eq__(const MolWallHitInfo& other) const;
  virtual bool eq_nonarray_attributes(const MolWallHitInfo& other, const bool ignore_name = false) const;
  bool operator == (const MolWallHitInfo& other) const { return __eq__(other);}
  bool operator != (const MolWallHitInfo& other) const { return !__eq__(other);}
  std::string to_str(const bool all_details=false, const std::string ind="") const ;

  // --- attributes ---
  int molecule_id;
  virtual void set_molecule_id(const int new_molecule_id_) {
    molecule_id = new_molecule_id_;
  }
  virtual int get_molecule_id() const {
    return molecule_id;
  }

  std::shared_ptr<GeometryObject> geometry_object;
  virtual void set_geometry_object(std::shared_ptr<GeometryObject> new_geometry_object_) {
    geometry_object = new_geometry_object_;
  }
  virtual std::shared_ptr<GeometryObject> get_geometry_object() const {
    return geometry_object;
  }

  int wall_index;
  virtual void set_wall_index(const int new_wall_index_) {
    wall_index = new_wall_index_;
  }
  virtual int get_wall_index() const {
    return wall_index;
  }

  double time;
  virtual void set_time(const double new_time_) {
    time = new_time_;
  }
  virtual double get_time() const {
    return time;
  }

  Vec3 pos3d;
  virtual void set_pos3d(const Vec3& new_pos3d_) {
    pos3d = new_pos3d_;
  }
  virtual const Vec3& get_pos3d() const {
    return pos3d;
  }

  double time_before_hit;
  virtual void set_time_before_hit(const double new_time_before_hit_) {
    time_before_hit = new_time_before_hit_;
  }
  virtual double get_time_before_hit() const {
    return time_before_hit;
  }

  Vec3 pos3d_before_hit;
  virtual void set_pos3d_before_hit(const Vec3& new_pos3d_before_hit_) {
    pos3d_before_hit = new_pos3d_before_hit_;
  }
  virtual const Vec3& get_pos3d_before_hit() const {
    return pos3d_before_hit;
  }

  // --- methods ---
}; // GenMolWallHitInfo

class MolWallHitInfo;
py::class_<MolWallHitInfo> define_pybinding_MolWallHitInfo(py::module& m);
} // namespace API
} // namespace MCell

#endif // API_GEN_MOL_WALL_HIT_INFO_H
