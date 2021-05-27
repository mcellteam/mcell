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

#ifndef API_GEN_MOLECULE_RELEASE_INFO_H
#define API_GEN_MOLECULE_RELEASE_INFO_H

#include "api/api_common.h"
#include "api/base_data_class.h"

namespace MCell {
namespace API {

class MoleculeReleaseInfo;
class Complex;
class PythonExportContext;

#define MOLECULE_RELEASE_INFO_CTOR() \
    MoleculeReleaseInfo( \
        std::shared_ptr<Complex> complex_, \
        const std::vector<double> location_ \
    ) { \
      class_name = "MoleculeReleaseInfo"; \
      complex = complex_; \
      location = location_; \
      postprocess_in_ctor(); \
      check_semantics(); \
    } \
    MoleculeReleaseInfo(DefaultCtorArgType) : \
      GenMoleculeReleaseInfo(DefaultCtorArgType()) { \
      set_all_attributes_as_default_or_unset(); \
    }

class GenMoleculeReleaseInfo: public BaseDataClass {
public:
  GenMoleculeReleaseInfo() {
  }
  GenMoleculeReleaseInfo(DefaultCtorArgType) {
  }
  void postprocess_in_ctor() override {}
  void check_semantics() const override;
  void set_initialized() override;
  void set_all_attributes_as_default_or_unset() override;

  MoleculeReleaseInfo copy_molecule_release_info() const;
  virtual bool __eq__(const MoleculeReleaseInfo& other) const;
  virtual bool eq_nonarray_attributes(const MoleculeReleaseInfo& other, const bool ignore_name = false) const;
  bool operator == (const MoleculeReleaseInfo& other) const { return __eq__(other);}
  bool operator != (const MoleculeReleaseInfo& other) const { return !__eq__(other);}
  std::string to_str(const std::string ind="") const override;

  std::string export_to_python(std::ostream& out, PythonExportContext& ctx) override;
  virtual std::string export_vec_location(std::ostream& out, PythonExportContext& ctx, const std::string& parent_name);


  // --- attributes ---
  std::shared_ptr<Complex> complex;
  virtual void set_complex(std::shared_ptr<Complex> new_complex_) {
    if (initialized) {
      throw RuntimeError("Value 'complex' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    complex = new_complex_;
  }
  virtual std::shared_ptr<Complex> get_complex() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return complex;
  }

  std::vector<double> location;
  virtual void set_location(const std::vector<double> new_location_) {
    if (initialized) {
      throw RuntimeError("Value 'location' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    location = new_location_;
  }
  virtual std::vector<double>& get_location() {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return location;
  }

  // --- methods ---
}; // GenMoleculeReleaseInfo

class MoleculeReleaseInfo;
py::class_<MoleculeReleaseInfo> define_pybinding_MoleculeReleaseInfo(py::module& m);
} // namespace API
} // namespace MCell

#endif // API_GEN_MOLECULE_RELEASE_INFO_H
