/******************************************************************************
 *
 * Copyright (C) 2020 by
 * The Salk Institute for Biological Studies
 *
 * Use of this source code is governed by an MIT-style
 * license that can be found in the LICENSE file or at
 * https://opensource.org/licenses/MIT.
 *
******************************************************************************/

#ifndef API_SPECIES_H
#define API_SPECIES_H

#include "generated/gen_species.h"
#include "generated/gen_constants.h"
#include "api/api_common.h"
#include "api/api_utils.h"
#include "api/complex.h"
#include "elementary_molecule_type.h"

#include "bng/bng_defines.h"

namespace MCell {
namespace API {

class Species: public GenSpecies {
public:
  SPECIES_CTOR()

  Species(Complex& cplx_inst)
    : GenSpecies(cplx_inst.to_bngl_str()) {

    set_all_custom_attributes_to_default();
    set_all_attributes_as_default_or_unset();

    // set copied value, cannot use GenSpecies ctor because we need to reset the attributes
    // warning: must be updated if Complex class' attributes change
    name = cplx_inst.name;
    elementary_molecules = cplx_inst.elementary_molecules;
    // should not be really used in Species but copying it as well for consistency
    orientation = cplx_inst.orientation;
  }

  void set_name(const std::string& name_) override {
    BaseDataClass::set_name(name_);
    // rerun initialization because the name is parsed as a BNGL string
    elementary_molecules.clear();
    postprocess_in_ctor();
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
    set_all_custom_attributes_to_default();

    // not calling derived check semantics
    if (get_num_set(name, elementary_molecules) != 1) {
      throw ValueError(
          S("Exactly one of ") + NAME_NAME + " or " + NAME_ELEMENTARY_MOLECULES +
          " must be set for " + NAME_CLASS_SPECIES + ".");
    }

    // 1) simple species defined by name
    if (is_set(name) && !is_set(elementary_molecules) &&
        (is_set(diffusion_constant_2d) || is_set(diffusion_constant_3d))) {

      if (!is_simple_species(name)) {
        throw ValueError("Only simple species can be fully defined by setting name and diffusion constant. "
            "For complex species, it cannot be usually deduced what the diffusion constants of elementary molecule types should be.");
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
      // will be added on initialization to model's elementary_molecule_types
      auto mt = std::make_shared<ElementaryMoleculeType>(
          name, std::vector<std::shared_ptr<ComponentType>>(),
          diffusion_constant_2d, diffusion_constant_3d,
          custom_time_step, custom_space_step,
          target_only
      );

      // and then molecule instance out of it
      elementary_molecules.push_back(
          std::make_shared<ElementaryMolecule>(mt)
      );
    }
    // 2) complex species defined through elementary_molecule_instances, or
    // 3) declaration (name is parsed in Complex::postprocess_in_ctor)
    else {
      // do semantic check
      check_no_extra_fields_are_set();
    }

    // need to finalize the initialization, also processes compartments
    Complex::postprocess_in_ctor();

    if (get_primary_compartment_name() == BNG::COMPARTMENT_NAME_IN ||
        get_primary_compartment_name() == BNG::COMPARTMENT_NAME_OUT) {
      throw ValueError(S(NAME_CLASS_SPECIES) + " with " + NAME_NAME + " " + name +
          " must not use compartments " + BNG::COMPARTMENT_NAME_IN + " or " + BNG::COMPARTMENT_NAME_OUT + ".");
    }
  }

  void set_all_custom_attributes_to_default() override {
    Complex::set_all_custom_attributes_to_default();
    species_id = SPECIES_ID_INVALID;
  }

  // using shorter printout when all_details is false
  std::string to_str(const bool all_details=false, const std::string ind="") const override;

  bool __eq__(const Species& other) const override;

  std::shared_ptr<Complex> inst(const Orientation orientation = Orientation::DEFAULT, const std::string& compartment_name = "") override {
    if (orientation != Orientation::DEFAULT && is_set(compartment_name)) {
      throw ValueError(S("Maximum one of ") + NAME_ORIENTATION + " or " + NAME_COMPARTMENT_NAME + " can be set not both.");
    }

    // make a deep copy and set extra attributes
    std::shared_ptr<Complex> res = deepcopy_complex();
    res->orientation = orientation;
    res->set_compartment_name(compartment_name);
    return res;
  }

  bool skip_python_export() const override {
    return is_species_superclass();
  }

  // default generated variant exports all details, we do not want this because
  // we export just the BNGL name and recompute the diffusion constant
  // and other fields from the elementary molecule types this species uses
  std::string export_to_python(std::ostream& out, PythonExportContext& ctx) override;

  bool is_species_superclass() const {
    return BNG::is_species_superclass(name);
  }

  bool warn_if_adding_identical_object() const { return true; }

  // simulation engine mapping
  species_id_t species_id;
};

} // namespace API
} // namespace MCell

#endif // API_SPECIES_H
