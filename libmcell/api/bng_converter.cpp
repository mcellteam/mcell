/******************************************************************************
 *
 * Copyright (C) 2020-2021 by
 * The Salk Institute for Biological Studies
 *
 * Use of this source code is governed by an MIT-style
 * license that can be found in the LICENSE file or at
 * https://opensource.org/licenses/MIT.
 *
******************************************************************************/

#include "api/bng_converter.h"

#include "api/component_type.h"
#include "api/elementary_molecule_type.h"
#include "api/complex.h"
#include "api/api_utils.h"

#include "bng/bng.h"

using namespace std;

namespace MCell {
namespace API {

BNG::component_type_id_t BNGConverter::convert_component_type(
    const std::string& elem_mol_type_name, API::ComponentType& api_ct) {
  // component types are identified by their name and set of allowed states, not just by their name
  BNG::ComponentType bng_ct;

  bng_ct.name = api_ct.name;
  bng_ct.elem_mol_type_name = elem_mol_type_name;

  for (string& s: api_ct.states) {
    bng_ct.allowed_state_ids.insert( bng_data.find_or_add_state_name(s) );
  }

  BNG::component_type_id_t ct_id = bng_data.find_or_add_component_type(bng_ct);
  if (ct_id == BNG::COMPONENT_TYPE_ID_INVALID) {
    throw ValueError("Conflicting component type " + api_ct.to_bngl_str() + ", " +
        " component with the same name already exists but has different allowed states.");
  }

  return ct_id;
}


BNG::Component BNGConverter::convert_component_instance(
    const std::string& elem_mol_type_name, API::Component& api_ci) {

  BNG::Component res(convert_component_type(elem_mol_type_name, *api_ci.component_type));

  if (api_ci.state == STATE_UNSET) {
    res.state_id = BNG::STATE_ID_DONT_CARE;
  }
  else {
    res.state_id = bng_data.find_state_id(api_ci.state);
    assert(res.state_id != BNG::STATE_ID_INVALID);
  }

  if (api_ci.bond == BOND_BOUND) {
    res.bond_value = BNG::BOND_VALUE_BOUND;
  }
  else if (api_ci.bond == BOND_UNBOUND) {
    res.bond_value = BNG::BOND_VALUE_UNBOUND;
  }
  else if (api_ci.bond == BOND_ANY) {
    res.bond_value = BNG::BOND_VALUE_ANY;
  }
  else {
    res.bond_value = api_ci.bond;
    assert(res.bond_value != BNG::BOND_VALUE_INVALID);
  }

  return res;
}


BNG::ElemMol BNGConverter::convert_molecule_instance(API::ElementaryMolecule& mi, const bool in_rxn_or_observables) {
  BNG::ElemMol res;

  res.elem_mol_type_id = convert_elementary_molecule_type(*mi.elementary_molecule_type, in_rxn_or_observables);

  for (std::shared_ptr<API::Component>& api_ci: mi.components) {
    res.components.push_back(convert_component_instance(mi.elementary_molecule_type->name, *api_ci));
  }

  if (is_set(mi.compartment_name)) {
    BNG::compartment_id_t comp_id = bng_data.find_compartment_id(mi.compartment_name);
    if (comp_id == BNG::COMPARTMENT_ID_INVALID) {
      throw ValueError("Elementary molecule " + mi.to_bngl_str() + " uses unknown compartment " + mi.compartment_name + ".");
    }
    res.compartment_id = comp_id;
  }
  else {
    res.compartment_id = BNG::COMPARTMENT_ID_NONE;
  }

  // we must also copy flags from the mol type
  res.finalize_flags_and_sort_components(bng_data);

  return res;
}


BNG::elem_mol_type_id_t BNGConverter::convert_elementary_molecule_type(
    API::ElementaryMoleculeType& api_mt, const bool in_rxn_or_observables) {
  if (api_mt.mol_type_id != BNG::ELEM_MOL_TYPE_ID_INVALID) {
    // already converted
    return api_mt.mol_type_id;
  }

  BNG::ElemMolType bng_mt;

  bng_mt.name = api_mt.name;

  if (!in_rxn_or_observables) {
    if (is_set(api_mt.diffusion_constant_2d)) {
      bng_mt.set_is_surf();
      bng_mt.D = api_mt.diffusion_constant_2d;
    }
    else if (is_set(api_mt.diffusion_constant_3d)) {
      bng_mt.set_is_vol();
      bng_mt.D = api_mt.diffusion_constant_3d;
    }
    else {
      throw RuntimeError(S("Diffusion constant for ") + NAME_CLASS_ELEMENTARY_MOLECULE_TYPE +
          " '" + bng_mt.name + "' was not set.");
    }

    if (is_set(api_mt.custom_time_step)) {
      bng_mt.custom_time_step = api_mt.custom_time_step;
    }
    else if (is_set(api_mt.custom_space_step)) {
      bng_mt.custom_space_step = api_mt.custom_space_step;
    }

    bng_mt.set_flag(BNG::SPECIES_MOL_FLAG_TARGET_ONLY, api_mt.target_only);
  }

  // components
  for (std::shared_ptr<API::ComponentType> api_ct: api_mt.components) {
    BNG::ComponentType bng_ct;

    bng_ct.name = api_ct->name;
    bng_ct.elem_mol_type_name = bng_mt.name;

    for (const string& state: api_ct->states) {
      bng_ct.allowed_state_ids.insert_unique(bng_data.find_or_add_state_name(state));
    }

    BNG::component_type_id_t ct_id = bng_data.find_or_add_component_type(bng_ct);
    if (ct_id == BNG::COMPONENT_TYPE_ID_INVALID) {
      throw ValueError("Conflicting component type " + api_ct->to_bngl_str() + ", " +
          " component with the same name already exists but has different allowed states.");
    }
    bng_mt.component_type_ids.push_back(ct_id);
  }

  if (!in_rxn_or_observables) {
    // we can convert only definitions
    bng_mt.compute_space_and_time_step(bng_config);
  }

  return bng_data.find_or_add_elem_mol_type(bng_mt);
}


BNG::Cplx BNGConverter::convert_complex(API::Complex& api_cplx, const bool in_observables, const bool in_rxn) {
  // create a temporary cplx instance that we will use for search
  BNG::Cplx bng_cplx(&bng_data);

  if (is_set(api_cplx.elementary_molecules)) {
    for (std::shared_ptr<API::ElementaryMolecule>& m: api_cplx.elementary_molecules) {
      BNG::ElemMol mi = convert_molecule_instance(*m, in_observables || in_rxn);

      bng_cplx.elem_mols.push_back(mi);
    }
  }
  else {
    throw ValueError(
        "Complex with name " + api_cplx.name + " does not have its " + NAME_ELEMENTARY_MOLECULES + " set. "
        "It should be always set because initialization of " + NAME_CLASS_COMPLEX + " from " + NAME_NAME + " is done "
        "in this class' constructor."
    );
  }

  // orientation or compartment does not have to be set for finalization,
  // this sets whether this is a surf or vol cplx
  bng_cplx.finalize_cplx();

  // BNG compartments were already created, they were also set for individual molecules
  if (is_set(api_cplx.compartment_name)) {
    // override all used compartments that were not set (are NONE)
    BNG::compartment_id_t in_out_id = BNG::get_in_or_out_compartment_id(api_cplx.compartment_name);
    if (in_out_id != BNG::COMPARTMENT_ID_INVALID) {
      bng_cplx.set_compartment_id(in_out_id, true);
    }
    else {
      const BNG::Compartment* bng_comp = bng_data.find_compartment(api_cplx.compartment_name);
      if (bng_cplx.is_vol() && (bng_comp == nullptr || !bng_comp->is_3d)) {
        throw ValueError("Did not find volume compartment " + api_cplx.compartment_name +
            " for a volume complex " + bng_cplx.to_str() + ".");
      }

      if (bng_cplx.is_surf() && (bng_comp == nullptr || bng_comp->is_3d)) {
        throw ValueError("Did not find surface compartment " + api_cplx.compartment_name +
            " for a surface complex " + bng_cplx.to_str() + ".");
      }

      bng_cplx.set_compartment_id(bng_comp->id, true);
    }
  }
  else {
    // main compartment was not set, do not change them

    if (!in_rxn && bng_cplx.is_vol() && api_cplx.orientation != Orientation::NONE && api_cplx.orientation != Orientation::DEFAULT) {
      throw ValueError("Orientation for a volume complex " + bng_cplx.to_str() +
          " must be set either to " + NAME_ENUM_ORIENTATION + "." + NAME_EV_NONE + " or " +
          NAME_ENUM_ORIENTATION + "." + NAME_EV_DEFAULT + ".");
    }
    else if (bng_cplx.is_surf() && api_cplx.orientation == Orientation::NONE) {
      throw ValueError("Orientation for a surface complex " + bng_cplx.to_str() +
          " must be set to a value other than " +  NAME_ENUM_ORIENTATION + "." + NAME_EV_NONE +
          " when " + NAME_COMPARTMENT_NAME + " is not specified.");
    }

    orientation_t orient = convert_api_orientation(api_cplx.orientation, true, bng_cplx.is_vol());
    bng_cplx.set_orientation(orient);
  }

  return bng_cplx;
}

} // namespace API
} // namespace MCell

