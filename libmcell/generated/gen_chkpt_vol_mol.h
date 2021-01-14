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

#ifndef API_GEN_CHKPT_VOL_MOL_H
#define API_GEN_CHKPT_VOL_MOL_H

#include "api/common.h"
#include "api\base_chkpt_mol.h"


namespace MCell {
namespace API {

class ChkptVolMol;
class Species;
class PythonExportContext;

#define CHKPT_VOL_MOL_CTOR() \
    ChkptVolMol( \
        const Vec3& pos_, \
        const int id_, \
        std::shared_ptr<Species> species_, \
        const float_t diffusion_time_, \
        const float_t birthday_, \
        const int flags_, \
        const float_t unimol_rx_time_ = FLT_UNSET \
    )  : GenChkptVolMol(id_,species_,diffusion_time_,birthday_,flags_,unimol_rx_time_) { \
      class_name = "ChkptVolMol"; \
      pos = pos_; \
      id = id_; \
      species = species_; \
      diffusion_time = diffusion_time_; \
      birthday = birthday_; \
      flags = flags_; \
      unimol_rx_time = unimol_rx_time_; \
      postprocess_in_ctor();\
      check_semantics();\
    }

class GenChkptVolMol: public BaseChkptMol {
public:
  GenChkptVolMol( 
      const int id_, 
      std::shared_ptr<Species> species_, 
      const float_t diffusion_time_, 
      const float_t birthday_, 
      const int flags_, 
      const float_t unimol_rx_time_ = FLT_UNSET 
  )  : BaseChkptMol(id_,species_,diffusion_time_,birthday_,flags_,unimol_rx_time_)  {
  }
  void postprocess_in_ctor() override {}
  void check_semantics() const override;
  void set_initialized() override;
  void set_all_attributes_as_default_or_unset() override;

  virtual bool __eq__(const ChkptVolMol& other) const;
  virtual bool eq_nonarray_attributes(const ChkptVolMol& other, const bool ignore_name = false) const;
  bool operator == (const ChkptVolMol& other) const { return __eq__(other);}
  bool operator != (const ChkptVolMol& other) const { return !__eq__(other);}
  std::string to_str(const std::string ind="") const override;

  virtual std::string export_to_python(std::ostream& out, PythonExportContext& ctx);


  // --- attributes ---
  Vec3 pos;
  virtual void set_pos(const Vec3& new_pos_) {
    if (initialized) {
      throw RuntimeError("Value 'pos' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    pos = new_pos_;
  }
  virtual const Vec3& get_pos() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return pos;
  }

  // --- methods ---
}; // GenChkptVolMol

class ChkptVolMol;
py::class_<ChkptVolMol> define_pybinding_ChkptVolMol(py::module& m);
} // namespace API
} // namespace MCell

#endif // API_GEN_CHKPT_VOL_MOL_H
