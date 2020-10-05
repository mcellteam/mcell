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

#ifndef API_SPECIES_H
#define API_SPECIES_H

#include <Python.h>

#include "generated/gen_species.h"
#include "generated/gen_constants.h"
#include "api/common.h"
#include "complex_instance.h"
#include "elementary_molecule_type.h"

#include "include/datamodel_defines.h"

namespace MCell {
namespace API {

class Species: public GenSpecies {
public:
  SPECIES_CTOR()

  // ctor for special ALL_*MOLECULES species
  // overlaps with generated ctor, the only difference is the type of the argument which is
  // const std::string&
  Species(const char* name_)
    : GenSpecies(name_) {
    set_all_attributes_as_default_or_unset();
    name = name_;
    species_id = SPECIES_ID_INVALID;
  }

  Species(ComplexInstance& cplx_inst)
    : GenSpecies(cplx_inst.to_bngl_str()),
      species_id(SPECIES_ID_INVALID) {
    set_all_attributes_as_default_or_unset();

    // set copied value, cannot use GenSpecies ctor because we need to reset the attributes
    // warning: must be updated if ComplexInstance attributes change
    name = cplx_inst.name;
    elementary_molecule_instances = cplx_inst.elementary_molecule_instances;
    // should not be really used in Species but copying it as well for consistency
    orientation = cplx_inst.orientation;
  }

  void postprocess_in_ctor() override {
    // initialization
    species_id = SPECIES_ID_INVALID;

    // we can do semantic checks also in postprocess_in_ctor
    if (elementary_molecule_instances.empty()) {
      if (!is_set(diffusion_constant_2d) && !is_set(diffusion_constant_3d)) {
        throw ValueError("Field diffusion_constant_2d or diffusion_constant_3d must be set for simple species.");
      }
      if (is_set(diffusion_constant_2d) && is_set(diffusion_constant_3d)) {
        throw ValueError("Only one of fields diffusion_constant_2d or diffusion_constant_3d can be set for simple species.");
      }

      // create a single ElementaryMoleculeType
      auto mt = std::make_shared<ElementaryMoleculeType>(
          name, std::vector<std::shared_ptr<ComponentType>>(),
          diffusion_constant_2d, diffusion_constant_3d
      );

      // and then molecule instance out of it
      elementary_molecule_instances.push_back(
          std::make_shared<ElementaryMoleculeInstance>(mt)
      );
    }
    else {
      // do semantic check
      if (is_set(diffusion_constant_2d)) {
        throw ValueError(S("Field ") + NAME_DIFFUSION_CONSTANT_2D +
            " must not be set for complex species because it is derived from its elementary molecule types.");
      }
      if (is_set(diffusion_constant_3d)) {
        throw ValueError(S("Field ") + NAME_DIFFUSION_CONSTANT_3D +
            " must not be set for complex species because it is derived from its elementary molecule types.");
      }
    }

    if (get_num_set(custom_time_step, custom_space_step) > 1) {
      throw ValueError(S("Only one of ") + NAME_CUSTOM_TIME_STEP + " or " + NAME_CUSTOM_SPACE_STEP +
          " may be set at the same time.");
    }
  }

  ComplexInstance inst(const Orientation orientation) override {
    // simply downcast because the possible definition of an underlying
    // ComplexInstance was done in postprocess_in_ctor
    // TODO: store orientation (?)
    ComplexInstance res = *dynamic_cast<ComplexInstance*>(this);
    res.orientation = orientation;
    return res;
  }

  bool is_species_superclass() const {
    return MCell::is_species_superclass(name);
  }

  // simulation engine mapping
  species_id_t species_id;
};

} // namespace API
} // namespace MCell

#endif // API_SPECIES_H
