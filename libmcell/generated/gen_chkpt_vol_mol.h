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

#ifndef API_GEN_CHKPT_VOL_MOL_H
#define API_GEN_CHKPT_VOL_MOL_H

#include "api/api_common.h"
#include "api/base_chkpt_mol.h"


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
        const double diffusion_time_, \
        const double birthday_, \
        const int flags_, \
        const double unimol_rxn_time_ = FLT_UNSET \
    )  : GenChkptVolMol(id_,species_,diffusion_time_,birthday_,flags_,unimol_rxn_time_) { \
      class_name = "ChkptVolMol"; \
      pos = pos_; \
      id = id_; \
      species = species_; \
      diffusion_time = diffusion_time_; \
      birthday = birthday_; \
      flags = flags_; \
      unimol_rxn_time = unimol_rxn_time_; \
      postprocess_in_ctor(); \
      check_semantics(); \
    } \
    ChkptVolMol(DefaultCtorArgType) : \
      GenChkptVolMol(DefaultCtorArgType()) { \
      set_all_attributes_as_default_or_unset(); \
      set_all_custom_attributes_to_default(); \
    }

class GenChkptVolMol: public BaseChkptMol {
public:
  GenChkptVolMol( 
      const int id_, 
      std::shared_ptr<Species> species_, 
      const double diffusion_time_, 
      const double birthday_, 
      const int flags_, 
      const double unimol_rxn_time_ = FLT_UNSET 
  )  : BaseChkptMol(id_,species_,diffusion_time_,birthday_,flags_,unimol_rxn_time_)  {
  }
  GenChkptVolMol() : BaseChkptMol(DefaultCtorArgType()) {
  }
  GenChkptVolMol(DefaultCtorArgType) {
  }
  void postprocess_in_ctor() override {}
  void check_semantics() const override;
  void set_initialized() override;
  void set_all_attributes_as_default_or_unset() override;

  std::shared_ptr<ChkptVolMol> copy_chkpt_vol_mol() const;
  std::shared_ptr<ChkptVolMol> deepcopy_chkpt_vol_mol(py::dict = py::dict()) const;
  virtual bool __eq__(const ChkptVolMol& other) const;
  virtual bool eq_nonarray_attributes(const ChkptVolMol& other, const bool ignore_name = false) const;
  bool operator == (const ChkptVolMol& other) const { return __eq__(other);}
  bool operator != (const ChkptVolMol& other) const { return !__eq__(other);}
  std::string to_str(const bool all_details=false, const std::string ind="") const override;

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
