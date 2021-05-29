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

#ifndef API_GEN_CHKPT_SURF_MOL_H
#define API_GEN_CHKPT_SURF_MOL_H

#include "api/api_common.h"
#include "api/base_chkpt_mol.h"


namespace MCell {
namespace API {

class ChkptSurfMol;
class GeometryObject;
class Species;
class PythonExportContext;

#define CHKPT_SURF_MOL_CTOR() \
    ChkptSurfMol( \
        const Vec2& pos_, \
        const Orientation orientation_, \
        std::shared_ptr<GeometryObject> geometry_object_, \
        const int wall_index_, \
        const int grid_tile_index_, \
        const int id_, \
        std::shared_ptr<Species> species_, \
        const double diffusion_time_, \
        const double birthday_, \
        const int flags_, \
        const double unimol_rx_time_ = FLT_UNSET \
    )  : GenChkptSurfMol(id_,species_,diffusion_time_,birthday_,flags_,unimol_rx_time_) { \
      class_name = "ChkptSurfMol"; \
      pos = pos_; \
      orientation = orientation_; \
      geometry_object = geometry_object_; \
      wall_index = wall_index_; \
      grid_tile_index = grid_tile_index_; \
      id = id_; \
      species = species_; \
      diffusion_time = diffusion_time_; \
      birthday = birthday_; \
      flags = flags_; \
      unimol_rx_time = unimol_rx_time_; \
      postprocess_in_ctor(); \
      check_semantics(); \
    } \
    ChkptSurfMol(DefaultCtorArgType) : \
      GenChkptSurfMol(DefaultCtorArgType()) { \
      set_all_attributes_as_default_or_unset(); \
      set_all_custom_attributes_to_default(); \
    }

class GenChkptSurfMol: public BaseChkptMol {
public:
  GenChkptSurfMol( 
      const int id_, 
      std::shared_ptr<Species> species_, 
      const double diffusion_time_, 
      const double birthday_, 
      const int flags_, 
      const double unimol_rx_time_ = FLT_UNSET 
  )  : BaseChkptMol(id_,species_,diffusion_time_,birthday_,flags_,unimol_rx_time_)  {
  }
  GenChkptSurfMol() : BaseChkptMol(DefaultCtorArgType()) {
  }
  GenChkptSurfMol(DefaultCtorArgType) {
  }
  void postprocess_in_ctor() override {}
  void check_semantics() const override;
  void set_initialized() override;
  void set_all_attributes_as_default_or_unset() override;

  std::shared_ptr<ChkptSurfMol> copy_chkpt_surf_mol() const;
  std::shared_ptr<ChkptSurfMol> deepcopy_chkpt_surf_mol(py::dict = py::dict()) const;
  virtual bool __eq__(const ChkptSurfMol& other) const;
  virtual bool eq_nonarray_attributes(const ChkptSurfMol& other, const bool ignore_name = false) const;
  bool operator == (const ChkptSurfMol& other) const { return __eq__(other);}
  bool operator != (const ChkptSurfMol& other) const { return !__eq__(other);}
  std::string to_str(const bool all_details=false, const std::string ind="") const override;

  virtual std::string export_to_python(std::ostream& out, PythonExportContext& ctx);


  // --- attributes ---
  Vec2 pos;
  virtual void set_pos(const Vec2& new_pos_) {
    if (initialized) {
      throw RuntimeError("Value 'pos' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    pos = new_pos_;
  }
  virtual const Vec2& get_pos() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return pos;
  }

  Orientation orientation;
  virtual void set_orientation(const Orientation new_orientation_) {
    if (initialized) {
      throw RuntimeError("Value 'orientation' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    orientation = new_orientation_;
  }
  virtual Orientation get_orientation() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return orientation;
  }

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

  int grid_tile_index;
  virtual void set_grid_tile_index(const int new_grid_tile_index_) {
    if (initialized) {
      throw RuntimeError("Value 'grid_tile_index' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    grid_tile_index = new_grid_tile_index_;
  }
  virtual int get_grid_tile_index() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return grid_tile_index;
  }

  // --- methods ---
}; // GenChkptSurfMol

class ChkptSurfMol;
py::class_<ChkptSurfMol> define_pybinding_ChkptSurfMol(py::module& m);
} // namespace API
} // namespace MCell

#endif // API_GEN_CHKPT_SURF_MOL_H
