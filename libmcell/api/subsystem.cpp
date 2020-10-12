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

#include "subsystem.h"

#include "bng/bng.h"

#include "api/component_type.h"
#include "api/component_instance.h"
#include "api/elementary_molecule_type.h"
#include "api/elementary_molecule_instance.h"
#include "api/complex.h"

#include "api/api_utils.h"

using namespace std;

namespace MCell {
namespace API {

const char* const MCELL_DIFFUSION_CONSTANT_3D_PREFIX = "MCELL_DIFFUSION_CONSTANT_3D_";
const char* const MCELL_DIFFUSION_CONSTANT_2D_PREFIX = "MCELL_DIFFUSION_CONSTANT_2D_";


void Subsystem::dump() const {
  std::cout << to_str() << "\n";
}


void Subsystem::load_bngl_molecule_types_and_reaction_rules(
    const std::string& file_name,
    const std::map<std::string, float_t>& parameter_overrides) {

  BNG::BNGData bng_data;

  int num_errors = BNG::parse_bngl_file(file_name, bng_data, parameter_overrides);
  if (num_errors != 0) {
    throw RuntimeError("Could not parse BNGL file " + file_name + ".");
  }

  // now convert everything we parsed into the API classes so that the user can
  // inspect or manipulate it if needed
  convert_bng_data_to_subsystem_data(bng_data);
}


void Subsystem::convert_bng_data_to_subsystem_data(const BNG::BNGData& bng_data) {
  // elementary molecules
  for (const BNG::MolType& mt: bng_data.get_molecule_types()) {
    convert_molecule_type(bng_data, mt);
  }

  // reaction rules
  for (const BNG::RxnRule& rr: bng_data.get_rxn_rules()) {
    convert_reaction_rule(bng_data, rr);
  }
}


void Subsystem::convert_molecule_type(const BNG::BNGData& bng_data, const BNG::MolType& bng_mt) {

  const string& name = bng_mt.name;
  auto res_mt = make_shared<API::ElementaryMoleculeType>(name);

  // using the MCELL_* parameters for now, see if the diffusion rate was
  // specified in the model
  float_t D2 = FLT_INVALID; // init to silence compiler warning
  bool found2 = bng_data.get_parameter_value(MCELL_DIFFUSION_CONSTANT_2D_PREFIX + name, D2);
  float_t D3 = FLT_INVALID;
  bool found3 = bng_data.get_parameter_value(MCELL_DIFFUSION_CONSTANT_3D_PREFIX + name, D3);

  if (found2 && found3) {
    throw RuntimeError("Molecule type '" + name + "' has both 2d and 3d diffusion constants specified.");
  }
  // no need to check whether the diffusion constant was set, it may be set later and
  // model conversion checks that it was set

  if (found2) {
    res_mt->diffusion_constant_2d = D2;
  }
  else if (found3) {
    res_mt->diffusion_constant_3d = D3;
  }

  // process components
  for (BNG::component_type_id_t ct_id: bng_mt.component_type_ids) {
    const BNG::ComponentType bng_ct = bng_data.get_component_type(ct_id);

    auto ct = make_shared<API::ComponentType>(bng_ct.name);

    // and allowed states
    for (BNG::state_id_t state_id: bng_ct.allowed_state_ids) {
      ct->states.push_back(bng_data.get_state_name(state_id));
    }

    res_mt->components.push_back(ct);
  }

  append_to_vec(elementary_molecule_types, res_mt);
}


void Subsystem::convert_reaction_rule(const BNG::BNGData& bng_data, const BNG::RxnRule& bng_rr) {

  // always a standard rxn
  if (bng_rr.type != BNG::RxnType::Standard) {
    throw RuntimeError("Unexpected type of reaction from BNGL file for reaction " +  bng_rr.name + ".");
  }

  auto res_rr = make_shared<API::ReactionRule>();
  if (bng_rr.name != "") {
    res_rr->name = bng_rr.name;
  }

  // RxnRule
  res_rr->fwd_rate = bng_rr.base_rate_constant;

  for (const BNG::Cplx& inst: bng_rr.reactants) {
    // MCell3R accepts reactants with any orientation
    res_rr->reactants.push_back(
        convert_cplx_instance_w_orientation(bng_data, inst, Orientation::ANY));
  }

  for (const BNG::Cplx& inst: bng_rr.products) {
    // MCell3R always creates products with the orientation up
    res_rr->products.push_back(
        convert_cplx_instance_w_orientation(bng_data, inst, Orientation::UP));
  }

  // allow reactions with identical names
  append_to_vec(reaction_rules, res_rr, false, true);
}


static int convert_bond_value(const BNG::bond_value_t bng_bond_value) {
  // we would like the BNG library to be independent, therefore
  // the BNG library does not use the same constants as the API
  switch (bng_bond_value) {
    case BNG::BOND_VALUE_INVALID:
      release_assert("Invalid bond value"); // should not really happen
      return -2;
    case BNG::BOND_VALUE_BOUND: // represents !+
      return BOND_BOUND;
    case BNG::BOND_VALUE_ANY: // represents !+
      return BOND_ANY;
    case BNG::BOND_VALUE_UNBOUND: // no
      return BOND_UNBOUND;
    default:
      return bng_bond_value;
  }
}



std::shared_ptr<API::Complex> Subsystem::convert_cplx_instance(
    const BNG::BNGData& bng_data,
    const BNG::Cplx& bng_cplx) {

  auto res_cplx_inst = API::Complex::make_shared_empty();

  // convert each molecule instance
  for (const BNG::MolInstance& bmg_mi: bng_cplx.mol_instances) {

    // find molecule type and create an instance
    const string& mt_name = bng_data.get_molecule_type(bmg_mi.mol_type_id).name;
    std::shared_ptr<ElementaryMoleculeType> api_emt = find_elementary_molecule_type(mt_name);
    assert(is_set(api_emt));

    // prepare a vector of component instances with their bonds set
    std::vector<std::shared_ptr<API::ComponentInstance>> api_comp_instances;
    for (const BNG::ComponentInstance& bng_ci: bmg_mi.component_instances) {
      const std::string& ct_name = bng_data.get_component_type(bng_ci.component_type_id).name;

      // we need to define component type here, they are not global
      // and belong to the elementary molecule type
      std::shared_ptr<API::ComponentType> api_ct = vec_find_by_name(api_emt->components, ct_name);
      auto api_comp_inst = make_shared<API::ComponentInstance>(api_ct);

      if (bng_ci.state_id != BNG::STATE_ID_DONT_CARE) {
        api_comp_inst->state = bng_data.get_state_name(bng_ci.state_id);
      }
      api_comp_inst->bond = convert_bond_value(bng_ci.bond_value);

      api_comp_instances.push_back(api_comp_inst);
    }

    // and append instantiated elementary molecule type
    res_cplx_inst->elementary_molecule_instances.push_back(api_emt->inst(api_comp_instances));
  }

  // set compartment
  if (bng_cplx.has_compartment()) {
    res_cplx_inst->compartment_name = bng_data.get_compartment(bng_cplx.get_compartment_id()).name;
  }

  return res_cplx_inst;
}


// sets orientation if the resulting cplx is a surface cplx
std::shared_ptr<API::Complex> Subsystem::convert_cplx_instance_w_orientation(
    const BNG::BNGData& bng_data,
    const BNG::Cplx& bng_inst,
    const Orientation orientation) {
  shared_ptr<API::Complex> res =
      Subsystem::convert_cplx_instance(bng_data, bng_inst);

  // only after conversion we can know whether a molecule is of surface volume type,
  // it is determined from the diffusion constants set with MCELL_DIFFUSION_CONSTANT_*
  if (res->is_surf()) {
    res->orientation = orientation;
  }

  return res;
}


} // namespace API
} // namespace MCell

