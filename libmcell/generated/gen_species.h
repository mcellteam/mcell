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

#ifndef API_GEN_SPECIES_H
#define API_GEN_SPECIES_H

#include "../api/common.h"
#include "../api/base_data_class.h"
#include "../api/complex_instance.h"

#include "../api/complex_instance.h"

namespace MCell {
namespace API {

class ComplexInstance;
class MoleculeInstance;

#define SPECIES_CTOR() \
    Species( \
        const std::string& name_, \
        const float_t diffusion_constant_2d_ = FLT_UNSET, \
        const float_t diffusion_constant_3d_ = FLT_UNSET, \
        const std::vector<std::shared_ptr<MoleculeInstance>> molecule_instances_ = std::vector<std::shared_ptr<MoleculeInstance>>(), \
        const Orientation orientation_ = Orientation::None \
    )  : GenSpecies(molecule_instances_,orientation_) { \
      class_name = "Species"; \
      name = name_; \
      diffusion_constant_2d = diffusion_constant_2d_; \
      diffusion_constant_3d = diffusion_constant_3d_; \
      molecule_instances = molecule_instances_; \
      orientation = orientation_; \
      postprocess_in_ctor();\
      check_semantics();\
    }

class GenSpecies: public ComplexInstance {
public:
  GenSpecies( 
      const std::vector<std::shared_ptr<MoleculeInstance>> molecule_instances_ = std::vector<std::shared_ptr<MoleculeInstance>>(), 
      const Orientation orientation_ = Orientation::None 
  )  : ComplexInstance(molecule_instances_,orientation_)  {
  }
  void postprocess_in_ctor() override {}
  void check_semantics() const override;
  void set_initialized() override;

  bool __eq__(const GenSpecies& other) const;
  std::string to_str(const std::string ind="") const override;

  // --- attributes ---
  float_t diffusion_constant_2d;
  virtual void set_diffusion_constant_2d(const float_t new_diffusion_constant_2d_) {
    if (initialized) {
      throw RuntimeError("Value 'diffusion_constant_2d' of object with name " + name + " (class " + class_name + ")"
                         "cannot be set after model was initialized.");
    }
    diffusion_constant_2d = new_diffusion_constant_2d_;
  }
  virtual float_t get_diffusion_constant_2d() const {
    return diffusion_constant_2d;
  }

  float_t diffusion_constant_3d;
  virtual void set_diffusion_constant_3d(const float_t new_diffusion_constant_3d_) {
    if (initialized) {
      throw RuntimeError("Value 'diffusion_constant_3d' of object with name " + name + " (class " + class_name + ")"
                         "cannot be set after model was initialized.");
    }
    diffusion_constant_3d = new_diffusion_constant_3d_;
  }
  virtual float_t get_diffusion_constant_3d() const {
    return diffusion_constant_3d;
  }

  // --- methods ---
  virtual ComplexInstance inst(const Orientation orientation = Orientation::NotSet) = 0;
}; // GenSpecies

class Species;
py::class_<Species> define_pybinding_Species(py::module& m);
} // namespace API
} // namespace MCell

#endif // API_GEN_SPECIES_H
