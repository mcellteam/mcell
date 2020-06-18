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

class Species;

#define MOLECULE_RELEASE_INFO_CTOR() \
    MoleculeReleaseInfo( \
        std::shared_ptr<Species> species_, \
        const std::vector<float_t> location_, \
        const Orientation orientation_ = Orientation::NONE \
    ) { \
      class_name = "MoleculeReleaseInfo"; \
      species = species_; \
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
  std::string to_str(const std::string ind="") const override;

  // --- attributes ---
  std::shared_ptr<Species> species;
  virtual void set_species(std::shared_ptr<Species> new_species_) {
    if (initialized) {
      throw RuntimeError("Value 'species' of object with name " + name + " (class " + class_name + ")"
                         "cannot be set after model was initialized.");
    }
    species = new_species_;
  }
  virtual std::shared_ptr<Species> get_species() const {
    return species;
  }

  std::vector<float_t> location;
  virtual void set_location(const std::vector<float_t> new_location_) {
    if (initialized) {
      throw RuntimeError("Value 'location' of object with name " + name + " (class " + class_name + ")"
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
      throw RuntimeError("Value 'orientation' of object with name " + name + " (class " + class_name + ")"
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
