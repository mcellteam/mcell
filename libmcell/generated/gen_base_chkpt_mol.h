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

#ifndef API_GEN_BASE_CHKPT_MOL_H
#define API_GEN_BASE_CHKPT_MOL_H

#include "api/api_common.h"
#include "api/base_data_class.h"

namespace MCell {
namespace API {

class BaseChkptMol;
class Species;
class PythonExportContext;

#define BASE_CHKPT_MOL_CTOR() \
    BaseChkptMol( \
        const int id_, \
        std::shared_ptr<Species> species_, \
        const double diffusion_time_, \
        const double birthday_, \
        const int flags_, \
        const double unimol_rxn_time_ = FLT_UNSET \
    ) { \
      class_name = "BaseChkptMol"; \
      id = id_; \
      species = species_; \
      diffusion_time = diffusion_time_; \
      birthday = birthday_; \
      flags = flags_; \
      unimol_rxn_time = unimol_rxn_time_; \
      postprocess_in_ctor(); \
      check_semantics(); \
    } \
    BaseChkptMol(DefaultCtorArgType) : \
      GenBaseChkptMol(DefaultCtorArgType()) { \
      set_all_attributes_as_default_or_unset(); \
      set_all_custom_attributes_to_default(); \
    }

class GenBaseChkptMol: public BaseDataClass {
public:
  GenBaseChkptMol() {
  }
  GenBaseChkptMol(DefaultCtorArgType) {
  }
  void postprocess_in_ctor() override {}
  void check_semantics() const override;
  void set_initialized() override;
  void set_all_attributes_as_default_or_unset() override;

  std::shared_ptr<BaseChkptMol> copy_base_chkpt_mol() const;
  std::shared_ptr<BaseChkptMol> deepcopy_base_chkpt_mol(py::dict = py::dict()) const;
  virtual bool __eq__(const BaseChkptMol& other) const;
  virtual bool eq_nonarray_attributes(const BaseChkptMol& other, const bool ignore_name = false) const;
  bool operator == (const BaseChkptMol& other) const { return __eq__(other);}
  bool operator != (const BaseChkptMol& other) const { return !__eq__(other);}
  std::string to_str(const bool all_details=false, const std::string ind="") const override;

  std::string export_to_python(std::ostream& out, PythonExportContext& ctx) override;


  // --- attributes ---
  int id;
  virtual void set_id(const int new_id_) {
    if (initialized) {
      throw RuntimeError("Value 'id' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    id = new_id_;
  }
  virtual int get_id() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return id;
  }

  std::shared_ptr<Species> species;
  virtual void set_species(std::shared_ptr<Species> new_species_) {
    if (initialized) {
      throw RuntimeError("Value 'species' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    species = new_species_;
  }
  virtual std::shared_ptr<Species> get_species() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return species;
  }

  double diffusion_time;
  virtual void set_diffusion_time(const double new_diffusion_time_) {
    if (initialized) {
      throw RuntimeError("Value 'diffusion_time' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    diffusion_time = new_diffusion_time_;
  }
  virtual double get_diffusion_time() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return diffusion_time;
  }

  double birthday;
  virtual void set_birthday(const double new_birthday_) {
    if (initialized) {
      throw RuntimeError("Value 'birthday' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    birthday = new_birthday_;
  }
  virtual double get_birthday() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return birthday;
  }

  int flags;
  virtual void set_flags(const int new_flags_) {
    if (initialized) {
      throw RuntimeError("Value 'flags' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    flags = new_flags_;
  }
  virtual int get_flags() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return flags;
  }

  double unimol_rxn_time;
  virtual void set_unimol_rxn_time(const double new_unimol_rxn_time_) {
    if (initialized) {
      throw RuntimeError("Value 'unimol_rxn_time' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    unimol_rxn_time = new_unimol_rxn_time_;
  }
  virtual double get_unimol_rxn_time() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return unimol_rxn_time;
  }

  // --- methods ---
}; // GenBaseChkptMol

class BaseChkptMol;
py::class_<BaseChkptMol> define_pybinding_BaseChkptMol(py::module& m);
} // namespace API
} // namespace MCell

#endif // API_GEN_BASE_CHKPT_MOL_H
