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

#ifndef API_GEN_RNG_STATE_H
#define API_GEN_RNG_STATE_H

#include "api/api_common.h"
#include "api/base_data_class.h"

namespace MCell {
namespace API {

class RngState;
class PythonExportContext;

#define RNG_STATE_CTOR() \
    RngState( \
        const uint64_t randcnt_, \
        const uint64_t aa_, \
        const uint64_t bb_, \
        const uint64_t cc_, \
        const std::vector<uint64_t> randslr_, \
        const std::vector<uint64_t> mm_, \
        const uint64_t rngblocks_ \
    ) { \
      class_name = "RngState"; \
      randcnt = randcnt_; \
      aa = aa_; \
      bb = bb_; \
      cc = cc_; \
      randslr = randslr_; \
      mm = mm_; \
      rngblocks = rngblocks_; \
      postprocess_in_ctor(); \
      check_semantics(); \
    } \
    RngState(DefaultCtorArgType) : \
      GenRngState(DefaultCtorArgType()) { \
      set_all_attributes_as_default_or_unset(); \
    }

class GenRngState: public BaseDataClass {
public:
  GenRngState() {
  }
  GenRngState(DefaultCtorArgType) {
  }
  void postprocess_in_ctor() override {}
  void check_semantics() const override;
  void set_initialized() override;
  void set_all_attributes_as_default_or_unset() override;

  std::shared_ptr<RngState> copy_rng_state() const;
  std::shared_ptr<RngState> deepcopy_rng_state(py::dict = py::dict()) const;
  virtual bool __eq__(const RngState& other) const;
  virtual bool eq_nonarray_attributes(const RngState& other, const bool ignore_name = false) const;
  bool operator == (const RngState& other) const { return __eq__(other);}
  bool operator != (const RngState& other) const { return !__eq__(other);}
  std::string to_str(const std::string ind="") const override;

  std::string export_to_python(std::ostream& out, PythonExportContext& ctx) override;
  virtual std::string export_vec_randslr(std::ostream& out, PythonExportContext& ctx, const std::string& parent_name);
  virtual std::string export_vec_mm(std::ostream& out, PythonExportContext& ctx, const std::string& parent_name);


  // --- attributes ---
  uint64_t randcnt;
  virtual void set_randcnt(const uint64_t new_randcnt_) {
    if (initialized) {
      throw RuntimeError("Value 'randcnt' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    randcnt = new_randcnt_;
  }
  virtual uint64_t get_randcnt() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return randcnt;
  }

  uint64_t aa;
  virtual void set_aa(const uint64_t new_aa_) {
    if (initialized) {
      throw RuntimeError("Value 'aa' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    aa = new_aa_;
  }
  virtual uint64_t get_aa() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return aa;
  }

  uint64_t bb;
  virtual void set_bb(const uint64_t new_bb_) {
    if (initialized) {
      throw RuntimeError("Value 'bb' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    bb = new_bb_;
  }
  virtual uint64_t get_bb() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return bb;
  }

  uint64_t cc;
  virtual void set_cc(const uint64_t new_cc_) {
    if (initialized) {
      throw RuntimeError("Value 'cc' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    cc = new_cc_;
  }
  virtual uint64_t get_cc() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return cc;
  }

  std::vector<uint64_t> randslr;
  virtual void set_randslr(const std::vector<uint64_t> new_randslr_) {
    if (initialized) {
      throw RuntimeError("Value 'randslr' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    randslr = new_randslr_;
  }
  virtual std::vector<uint64_t>& get_randslr() {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return randslr;
  }

  std::vector<uint64_t> mm;
  virtual void set_mm(const std::vector<uint64_t> new_mm_) {
    if (initialized) {
      throw RuntimeError("Value 'mm' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    mm = new_mm_;
  }
  virtual std::vector<uint64_t>& get_mm() {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return mm;
  }

  uint64_t rngblocks;
  virtual void set_rngblocks(const uint64_t new_rngblocks_) {
    if (initialized) {
      throw RuntimeError("Value 'rngblocks' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    rngblocks = new_rngblocks_;
  }
  virtual uint64_t get_rngblocks() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return rngblocks;
  }

  // --- methods ---
}; // GenRngState

class RngState;
py::class_<RngState> define_pybinding_RngState(py::module& m);
} // namespace API
} // namespace MCell

#endif // API_GEN_RNG_STATE_H
