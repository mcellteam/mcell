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

#ifndef API_SPECIES_H
#define API_SPECIES_H

#include <string>

#include "../generated/gen_species.h"
#include "common.h"
#include "complex_instance.h"
#include "molecule_type.h"

namespace MCell {
namespace API {

class Species: public GenSpecies {
public:
  // having a generated ctor makes changes much simpler
  SPECIES_CTOR()

  void ctor_postprocess() override {
    // something like a semantic check? do semantic check in constructor?
    if (molecule_types.empty()) {
      if (!is_set(diffusion_constant_2d) && !is_set(diffusion_constant_3d)) {
        throw SemanticException("Field diffusion_constant_2d or diffusion_constant_3d must be set for simple species.");
      }

      // create a single MoleculeType
      auto mt = std::make_shared<MoleculeType>(
          name, std::vector<std::shared_ptr<ComponentType>>(),
          diffusion_constant_2d, diffusion_constant_3d
      );

      // and then molecule instance out of it
      molecule_types.push_back(
          std::make_shared<MoleculeInstance>(mt)
      );
    }
    else {
      // do semantic check
      if (is_set(diffusion_constant_2d)) {
        throw SemanticException("Field diffusion_constant_2d must not be set for simple species.");
      }
      if (is_set(diffusion_constant_3d)) {
        throw SemanticException("Field diffusion_constant_3d must not be set for simple species.");
      }
    }
  }

  ComplexInstance inst(const int orientation) override {
    // simply downcast because the possible definition of an underlying
    // ComplexInstance was done in ctor_postprocess
    return *dynamic_cast<ComplexInstance*>(this);
  }

};

} // namespace API
} // namespace MCell

#endif // API_SPECIES_H
