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

#ifndef API_GEN_BASE_CHKPT_MOL_H
#define API_GEN_BASE_CHKPT_MOL_H

#include "api/common.h"
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
        const float_t diffusion_time_, \
        const float_t unimol_rx_time_, \
        const float_t birthday_ \
    ) { \
      class_name = "BaseChkptMol"; \
      id = id_; \
      species = species_; \
      diffusion_time = diffusion_time_; \
      unimol_rx_time = unimol_rx_time_; \
      birthday = birthday_; \
      postprocess_in_ctor();\
      check_semantics();\
    }

class GenBaseChkptMol: public BaseDataClass {
public:
  void postprocess_in_ctor() override {}
  void check_semantics() const override;
  void set_initialized() override;
  void set_all_attributes_as_default_or_unset() override;

  virtual bool __eq__(const BaseChkptMol& other) const;
  virtual bool eq_nonarray_attributes(const BaseChkptMol& other, const bool ignore_name = false) const;
  bool operator == (const BaseChkptMol& other) const { return __eq__(other);}
  bool operator != (const BaseChkptMol& other) const { return !__eq__(other);}
  std::string to_str(const std::string ind="") const override;

  std::string export_to_python(std::ostream& out, PythonExportContext& ctx) const override;


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

  float_t diffusion_time;
  virtual void set_diffusion_time(const float_t new_diffusion_time_) {
    if (initialized) {
      throw RuntimeError("Value 'diffusion_time' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    diffusion_time = new_diffusion_time_;
  }
  virtual float_t get_diffusion_time() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return diffusion_time;
  }

  float_t unimol_rx_time;
  virtual void set_unimol_rx_time(const float_t new_unimol_rx_time_) {
    if (initialized) {
      throw RuntimeError("Value 'unimol_rx_time' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    unimol_rx_time = new_unimol_rx_time_;
  }
  virtual float_t get_unimol_rx_time() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return unimol_rx_time;
  }

  float_t birthday;
  virtual void set_birthday(const float_t new_birthday_) {
    if (initialized) {
      throw RuntimeError("Value 'birthday' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    birthday = new_birthday_;
  }
  virtual float_t get_birthday() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return birthday;
  }

  // --- methods ---
}; // GenBaseChkptMol

class BaseChkptMol;
py::class_<BaseChkptMol> define_pybinding_BaseChkptMol(py::module& m);
} // namespace API
} // namespace MCell

#endif // API_GEN_BASE_CHKPT_MOL_H
