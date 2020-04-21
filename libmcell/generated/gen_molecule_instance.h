/******************************************************************************
 *
 * Copyright (C) 2020 by
 * The Salk Institute for Biological Studies and
 * Pittsburgh Supercomputing Center, Carnegie Mellon University
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

#ifndef API_GEN_MOLECULE_INSTANCE_H
#define API_GEN_MOLECULE_INSTANCE_H

#include "../api/common.h"

namespace MCell {
namespace API {

class ComponentInstance;
class MoleculeType;

#define MOLECULE_INSTANCE_CTOR() \
    MoleculeInstance( \
        const MoleculeType* molecule_type_, \
        const std::vector<ComponentInstance*> components_ = std::vector<ComponentInstance*>() \
    ) { \
      class_name = "MoleculeInstance"; \
      molecule_type = molecule_type_; \
      components = components_; \
    }

class GenMoleculeInstance: public BaseDataClass {
public:
  SemRes check_semantics(std::ostream& out) const override;
  std::string to_str() const override;

  // --- attributes ---
  const MoleculeType* molecule_type;
  virtual void set_molecule_type(const MoleculeType* new_molecule_type_) {
    molecule_type = new_molecule_type_;
  }
  virtual const MoleculeType* get_molecule_type() const {
    return molecule_type;
  }

  std::vector<ComponentInstance*> components;
  virtual void set_components(const std::vector<ComponentInstance*> new_components_) {
    components = new_components_;
  }
  virtual std::vector<ComponentInstance*> get_components() const {
    return components;
  }

  // --- methods ---
}; // GenMoleculeInstance

class MoleculeInstance;
py::class_<MoleculeInstance> define_pybinding_MoleculeInstance(py::module& m);
} // namespace API
} // namespace MCell

#endif // API_GEN_MOLECULE_INSTANCE_H
