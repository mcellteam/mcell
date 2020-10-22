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

#include "../api/common.h"
#include "../api/base_data_class.h"

namespace MCell {
namespace API {

class Complex;

#define MOLECULE_RELEASE_INFO_CTOR() \
    MoleculeReleaseInfo( \
        std::shared_ptr<Complex> complex_, \
        const std::vector<float_t> location_, \
        const Orientation orientation_ = Orientation::NONE \
    ) { \
      class_name = "MoleculeReleaseInfo"; \
      complex = complex_; \
      location = location_; \
      orientation = orientation_; \
      postprocess_in_ctor();\
      check_semantics();\
    }

class GenMoleculeReleaseInfo: public BaseDataClass {
public:
  void postprocess_in_ctor() override {}
  void check_semantics() const override;
  void set_initialized() override;
  bool __eq__(const GenMoleculeReleaseInfo& other) const;
  void set_all_attributes_as_default_or_unset() override;

  std::string to_str(const std::string ind="") const override;

  // --- attributes ---
  std::shared_ptr<Complex> complex;
  virtual void set_complex(std::shared_ptr<Complex> new_complex_) {
    if (initialized) {
      throw RuntimeError("Value 'complex' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    complex = new_complex_;
  }
  virtual std::shared_ptr<Complex> get_complex() const {
    return complex;
  }

  std::vector<float_t> location;
  virtual void set_location(const std::vector<float_t> new_location_) {
    if (initialized) {
      throw RuntimeError("Value 'location' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    location = new_location_;
  }
  virtual std::vector<float_t> get_location() const {
    return location;
  }

  Orientation orientation;
  virtual void set_orientation(const Orientation new_orientation_) {
    if (initialized) {
      throw RuntimeError("Value 'orientation' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    orientation = new_orientation_;
  }
  virtual Orientation get_orientation() const {
    return orientation;
  }

  // --- methods ---
}; // GenMoleculeReleaseInfo

class MoleculeReleaseInfo;
py::class_<MoleculeReleaseInfo> define_pybinding_MoleculeReleaseInfo(py::module& m);
} // namespace API
} // namespace MCell

#endif // API_GEN_MOLECULE_RELEASE_INFO_H
