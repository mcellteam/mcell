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
#include "api/api_utils.h"
#include "api/complex.h"
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

  Species(Complex& cplx_inst)
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

  void check_no_extra_fields_are_set() const {
    std::string msg =
        " must not be set for complex species because it is derived from its elementary molecule types.";
    if (is_set(diffusion_constant_2d)) {
      throw ValueError(S("Field ") + NAME_DIFFUSION_CONSTANT_2D + msg);
    }
    if (is_set(diffusion_constant_3d)) {
      throw ValueError(S("Field ") + NAME_DIFFUSION_CONSTANT_3D + msg);
    }
    if (is_set(custom_time_step)) {
      throw ValueError(S("Field ") + NAME_CUSTOM_TIME_STEP + msg);
    }
    if (is_set(custom_space_step)) {
      throw ValueError(S("Field ") + NAME_CUSTOM_SPACE_STEP + msg);
    }
    if (target_only) {
      throw ValueError(S("Field ") + NAME_TARGET_ONLY + msg);
    }
  }


  // we are making changes, so semantic checks are here instead of in const check_semantics
  void postprocess_in_ctor() override {
    // initialization
    species_id = SPECIES_ID_INVALID;

    // not calling derived check semantics
    if (get_num_set(name, elementary_molecule_instances) != 1) {
      throw ValueError(
          S("Exactly one of ") + NAME_NAME + " or " + NAME_ELEMENTARY_MOLECULE_INSTANCES +
          " must be set for " + NAME_CLASS_SPECIES + ".");
    }

    // 1) simple species defined by name
    if (is_set(name) && !is_set(elementary_molecule_instances) && (is_set(diffusion_constant_2d) || is_set(diffusion_constant_3d))) {

      if (!is_simple_species(name)) {
        throw ValueError("Only simple species can be fully defined by setting name and diffusion constant.");
      }

      if (name.find('.') != std::string::npos) {
        throw ValueError("Simple species name must not contain '.', this is incompatible with BNGL definition"
            ", error for " + name + ".");
      }

      if (is_set(diffusion_constant_2d) && is_set(diffusion_constant_3d)) {
        throw ValueError("Only one of fields diffusion_constant_2d or diffusion_constant_3d can be set for simple species.");
      }

      if (get_num_set(custom_time_step, custom_space_step) > 1) {
        throw ValueError(S("Only one of ") + NAME_CUSTOM_TIME_STEP + " or " + NAME_CUSTOM_SPACE_STEP +
            " may be set at the same time.");
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
    // 2) complex species defined through elementary_molecule_instances
    else if (!elementary_molecule_instances.empty()) {
      // do semantic check
      check_no_extra_fields_are_set();
    }
    // 3) declaration
    else {
      // TODO: we can check that the BNGL string is correct here
      check_no_extra_fields_are_set();
    }
  }

  Complex inst(const Orientation orientation = Orientation::DEFAULT, const std::string& compartment_name = "") override {
    if (orientation != Orientation::DEFAULT && is_set(compartment_name)) {
      throw ValueError(S("Maximum one of ") + NAME_ORIENTATION + " or " + NAME_COMPARTMENT_NAME + " can be set not both.");
    }

    // simply downcast to Complex and set extra attributes
    Complex res = *dynamic_cast<Complex*>(this);
    res.orientation = orientation;
    res.compartment_name = compartment_name;
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
