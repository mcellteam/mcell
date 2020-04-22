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

#ifndef API_GEN_COMPLEX_INSTANCE_H
#define API_GEN_COMPLEX_INSTANCE_H

#include "../api/common.h"

namespace MCell {
namespace API {

class MoleculeInstance;

#define COMPLEX_INSTANCE_CTOR() \
    ComplexInstance( \
        const std::vector<std::shared_ptr<MoleculeInstance>> molecule_types_ = std::vector<std::shared_ptr<MoleculeInstance>>() \
    ) { \
      class_name = "ComplexInstance"; \
      molecule_types = molecule_types_; \
    }

class GenComplexInstance: public BaseDataClass {
public:
  SemRes check_semantics(std::ostream& out) const override;
  std::string to_str(const std::string ind="") const override;

  // --- attributes ---
  std::vector<std::shared_ptr<MoleculeInstance>> molecule_types;
  virtual void set_molecule_types(const std::vector<std::shared_ptr<MoleculeInstance>> new_molecule_types_) {
    molecule_types = new_molecule_types_;
  }
  virtual std::vector<std::shared_ptr<MoleculeInstance>> get_molecule_types() const {
    return molecule_types;
  }

  // --- methods ---
}; // GenComplexInstance

class ComplexInstance;
py::class_<ComplexInstance> define_pybinding_ComplexInstance(py::module& m);
} // namespace API
} // namespace MCell

#endif // API_GEN_COMPLEX_INSTANCE_H