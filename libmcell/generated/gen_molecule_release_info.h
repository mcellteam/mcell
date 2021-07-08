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
      set_all_custom_attributes_to_default(); \
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

  std::shared_ptr<MoleculeReleaseInfo> copy_molecule_release_info() const;
  std::shared_ptr<MoleculeReleaseInfo> deepcopy_molecule_release_info(py::dict = py::dict()) const;
  virtual bool __eq__(const MoleculeReleaseInfo& other) const;
  virtual bool eq_nonarray_attributes(const MoleculeReleaseInfo& other, const bool ignore_name = false) const;
  bool operator == (const MoleculeReleaseInfo& other) const { return __eq__(other);}
  bool operator != (const MoleculeReleaseInfo& other) const { return !__eq__(other);}
  std::string to_str(const bool all_details=false, const std::string ind="") const override;

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
