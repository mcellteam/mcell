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

#ifndef API_GEN_COMPLEX_INSTANCE_H
#define API_GEN_COMPLEX_INSTANCE_H

#include "../api/common.h"
#include "../api/base_data_class.h"

namespace MCell {
namespace API {

class MoleculeInstance;

#define COMPLEX_INSTANCE_CTOR() \
    ComplexInstance( \
        const std::vector<std::shared_ptr<MoleculeInstance>> molecule_instances_ = std::vector<std::shared_ptr<MoleculeInstance>>(), \
        const Orientation orientation_ = Orientation::None \
    ) { \
      class_name = "ComplexInstance"; \
      molecule_instances = molecule_instances_; \
      orientation = orientation_; \
      postprocess_in_ctor();\
      check_semantics();\
    }

class GenComplexInstance: public BaseDataClass {
public:
  void postprocess_in_ctor() override {}
  void check_semantics() const override;
  void set_initialized() override;

  bool __eq__(const GenComplexInstance& other) const;
  std::string to_str(const std::string ind="") const override;

  // --- attributes ---
  std::vector<std::shared_ptr<MoleculeInstance>> molecule_instances;
  virtual void set_molecule_instances(const std::vector<std::shared_ptr<MoleculeInstance>> new_molecule_instances_) {
    if (initialized) {
      throw RuntimeError("Value 'molecule_instances' of object with name " + name + " (class " + class_name + ")"
                         "cannot be set after model was initialized.");
    }
    molecule_instances = new_molecule_instances_;
  }
  virtual std::vector<std::shared_ptr<MoleculeInstance>> get_molecule_instances() const {
    return molecule_instances;
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
}; // GenComplexInstance

class ComplexInstance;
py::class_<ComplexInstance> define_pybinding_ComplexInstance(py::module& m);
} // namespace API
} // namespace MCell

#endif // API_GEN_COMPLEX_INSTANCE_H