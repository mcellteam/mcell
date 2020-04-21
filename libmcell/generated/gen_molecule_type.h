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

#ifndef API_GEN_MOLECULE_TYPE_H
#define API_GEN_MOLECULE_TYPE_H

#include "../api/common.h"

namespace MCell {
namespace API {

class ComponentInstance;
class ComponentType;
class MoleculeInstance;

#define MOLECULE_TYPE_CTOR() \
    MoleculeType( \
        const std::string& name_, \
        const std::vector<std::shared_ptr<ComponentType>> components_ = std::vector<std::shared_ptr<ComponentType>>(), \
        const float_t diffusion_constant_2d_ = FLT_UNSET, \
        const float_t diffusion_constant_3d_ = FLT_UNSET \
    ) { \
      class_name = "MoleculeType"; \
      name = name_; \
      components = components_; \
      diffusion_constant_2d = diffusion_constant_2d_; \
      diffusion_constant_3d = diffusion_constant_3d_; \
    }

class GenMoleculeType: public BaseDataClass {
public:
  SemRes check_semantics(std::ostream& out) const override;
  std::string to_str() const override;

  // --- attributes ---
  std::vector<std::shared_ptr<ComponentType>> components;
  virtual void set_components(const std::vector<std::shared_ptr<ComponentType>> new_components_) {
    components = new_components_;
  }
  virtual std::vector<std::shared_ptr<ComponentType>> get_components() const {
    return components;
  }

  float_t diffusion_constant_2d;
  virtual void set_diffusion_constant_2d(const float_t new_diffusion_constant_2d_) {
    diffusion_constant_2d = new_diffusion_constant_2d_;
  }
  virtual float_t get_diffusion_constant_2d() const {
    return diffusion_constant_2d;
  }

  float_t diffusion_constant_3d;
  virtual void set_diffusion_constant_3d(const float_t new_diffusion_constant_3d_) {
    diffusion_constant_3d = new_diffusion_constant_3d_;
  }
  virtual float_t get_diffusion_constant_3d() const {
    return diffusion_constant_3d;
  }

  // --- methods ---
  virtual MoleculeInstance inst(const std::vector<std::shared_ptr<ComponentInstance>> components) = 0;
}; // GenMoleculeType

class MoleculeType;
py::class_<MoleculeType> define_pybinding_MoleculeType(py::module& m);
} // namespace API
} // namespace MCell

#endif // API_GEN_MOLECULE_TYPE_H
