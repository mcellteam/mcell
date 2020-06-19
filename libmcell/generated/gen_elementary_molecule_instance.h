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

#ifndef API_GEN_ELEMENTARY_MOLECULE_INSTANCE_H
#define API_GEN_ELEMENTARY_MOLECULE_INSTANCE_H

#include "../api/common.h"
#include "../api/base_data_class.h"

namespace MCell {
namespace API {

class ComponentInstance;
class ElementaryMoleculeType;

#define ELEMENTARY_MOLECULE_INSTANCE_CTOR() \
    ElementaryMoleculeInstance( \
        std::shared_ptr<ElementaryMoleculeType> elementary_molecule_type_, \
        const std::vector<std::shared_ptr<ComponentInstance>> components_ = std::vector<std::shared_ptr<ComponentInstance>>() \
    ) { \
      class_name = "ElementaryMoleculeInstance"; \
      elementary_molecule_type = elementary_molecule_type_; \
      components = components_; \
      postprocess_in_ctor();\
      check_semantics();\
    }

class GenElementaryMoleculeInstance: public BaseDataClass {
public:
  void postprocess_in_ctor() override {}
  void check_semantics() const override;
  void set_initialized() override;

  bool __eq__(const GenElementaryMoleculeInstance& other) const;
  std::string to_str(const std::string ind="") const override;

  // --- attributes ---
  std::shared_ptr<ElementaryMoleculeType> elementary_molecule_type;
  virtual void set_elementary_molecule_type(std::shared_ptr<ElementaryMoleculeType> new_elementary_molecule_type_) {
    if (initialized) {
      throw RuntimeError("Value 'elementary_molecule_type' of object with name " + name + " (class " + class_name + ")"
                         "cannot be set after model was initialized.");
    }
    elementary_molecule_type = new_elementary_molecule_type_;
  }
  virtual std::shared_ptr<ElementaryMoleculeType> get_elementary_molecule_type() const {
    return elementary_molecule_type;
  }

  std::vector<std::shared_ptr<ComponentInstance>> components;
  virtual void set_components(const std::vector<std::shared_ptr<ComponentInstance>> new_components_) {
    if (initialized) {
      throw RuntimeError("Value 'components' of object with name " + name + " (class " + class_name + ")"
                         "cannot be set after model was initialized.");
    }
    components = new_components_;
  }
  virtual std::vector<std::shared_ptr<ComponentInstance>> get_components() const {
    return components;
  }

  // --- methods ---
  virtual std::string to_bngl_str() = 0;
}; // GenElementaryMoleculeInstance

class ElementaryMoleculeInstance;
py::class_<ElementaryMoleculeInstance> define_pybinding_ElementaryMoleculeInstance(py::module& m);
} // namespace API
} // namespace MCell

#endif // API_GEN_ELEMENTARY_MOLECULE_INSTANCE_H
