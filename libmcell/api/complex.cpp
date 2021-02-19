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

#include "api/complex.h"
#include "api/species.h"
#include "api/subsystem.h"
#include "api/component.h"
#include "api/component_type.h"
#include "api/elementary_molecule.h"
#include "api/elementary_molecule_type.h"
#include "bng/bng.h"

using namespace std;

namespace MCell {
namespace API {



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


// make a deep copy, used from Species::inst
std::shared_ptr<Complex> Complex::clone() const {
  std::shared_ptr<Complex> res = make_shared<Complex>();

  res->name = name;
  res->orientation = orientation;
  res->compartment_name = compartment_name;

  for (const auto& em: elementary_molecules) {
    res->elementary_molecules.push_back(em->clone());
  }

  return res;
}


std::shared_ptr<API::Complex> Complex::construct_from_bng_cplx(
    const BNG::BNGData& bng_data,
    const BNG::Cplx& bng_cplx) {

  std::shared_ptr<API::Complex> res_cplx_inst = API::Complex::construct_empty();

  // convert each molecule instance
  for (const BNG::ElemMol& bmg_em: bng_cplx.elem_mols) {

    // find molecule type and create an instance
    const BNG::ElemMolType& emt = bng_data.get_elem_mol_type(bmg_em.elem_mol_type_id);
    std::shared_ptr<ElementaryMoleculeType> api_emt =
        ElementaryMoleculeType::construct_from_bng_elem_mol_type(bng_data, emt);
    assert(is_set(api_emt));

    // prepare a vector of component instances with their bonds set
    std::vector<std::shared_ptr<API::Component>> api_comp_instances;
    for (const BNG::Component& bng_ci: bmg_em.components) {
      const std::string& ct_name = bng_data.get_component_type(bng_ci.component_type_id).name;

      // we need to define component type here, they are not global
      // and belong to the elementary molecule type
      std::shared_ptr<API::ComponentType> api_ct = vec_find_by_name(api_emt->components, ct_name);
      auto api_comp_inst = make_shared<API::Component>(api_ct);

      if (bng_ci.state_id != BNG::STATE_ID_DONT_CARE) {
        api_comp_inst->state = bng_data.get_state_name(bng_ci.state_id);
      }
      api_comp_inst->bond = convert_bond_value(bng_ci.bond_value);

      api_comp_instances.push_back(api_comp_inst);
    }

    // determine compartment name
    BNG::compartment_id_t cid = bmg_em.compartment_id;
    string compartment_name;
    if (cid == BNG::COMPARTMENT_ID_NONE) {
      compartment_name = STR_UNSET;
    }
    else if (BNG::is_in_out_compartment_id(cid)) {
      compartment_name = BNG::compartment_id_to_str(cid);
    }
    else {
      compartment_name = bng_data.get_compartment(cid).name;
    }

    // and append instantiated elementary molecule type
    res_cplx_inst->elementary_molecules.push_back(api_emt->inst(api_comp_instances, compartment_name));
  }

  return res_cplx_inst;
}


// sets orientation if the resulting cplx is a surface cplx
std::shared_ptr<API::Complex> Complex::construct_from_bng_cplx_w_orientation(
    const BNG::BNGData& bng_data,
    const BNG::Cplx& bng_inst,
    const Orientation orientation) {
  shared_ptr<API::Complex> res =
      Complex::construct_from_bng_cplx(bng_data, bng_inst);

  // only after conversion we can know whether a molecule is of surface volume type,
  // it is determined from the diffusion constants set with MCELL_DIFFUSION_CONSTANT_*
  if (res->is_surf()) {
    res->orientation = orientation;
  }

  return res;
}


void Complex::postprocess_in_ctor() {
  if (name == STR_UNSET) {
    // ignore
    return;
  }

  // species did its own pre-initialization and simple complex has its
  // elementary molecule fully defined
  if (!is_species_object() && get_num_set(name, elementary_molecules) != 1) {
    throw ValueError(
        S("Exactly one of ") + NAME_NAME + " or " + NAME_ELEMENTARY_MOLECULES +
        " must be set for " + NAME_CLASS_COMPLEX + ".");
  }

  if (is_set(name)) {
    if (is_simple_species(name) && name.find('.') != std::string::npos) {
      throw ValueError("Simple species or complex name must not contain '.', this is incompatible with BNGL definition"
          ", error for " + name + ".");
    }
  }

  // if name was set, parse it and initialize elementary_molecules
  if (is_set(name) && !is_set(elementary_molecules)) {
    // parse BNGL string
    BNG::BNGData local_bng_data;
    BNG::Cplx bng_cplx(&local_bng_data);
    int num_errors = BNG::parse_single_cplx_string(name, local_bng_data, bng_cplx);
    if (num_errors) {
      throw ValueError("Could not parse BNGL string " + name + " that defines a " + NAME_CLASS_COMPLEX + ".");
    }

    // convert BNG data to elementary_molecules
    std::shared_ptr<API::Complex> converted_complex =
        Complex::construct_from_bng_cplx(local_bng_data, bng_cplx);

    // set compartment name - must be done before elementary molecules are added
    if (is_set(converted_complex->compartment_name)) {
      if (is_set(compartment_name) && compartment_name != converted_complex->compartment_name) {
        throw ValueError("Inconsistent compartment usage in complex " + name +
            ", parsed compartment is " + converted_complex->compartment_name + " but compartment set as " +
            NAME_COMPARTMENT_NAME + " is " + compartment_name + ".");
      }
      set_compartment_name(converted_complex->compartment_name);
    }

    // copy resulting elementary molecules array - their compartments are set
    elementary_molecules = converted_complex->elementary_molecules;
  }

  check_and_set_compartments();

  assert(!elementary_molecules.empty());
}


const std::string& Complex::get_canonical_name() const {
  if (cached_data_are_uptodate && canonical_name != "") {
    return canonical_name;
  }

  // note: using the elementary_molecules array would be faster but using parser is simpler to implement

  // get BNGL string and parse it, then get canonical name
  BNG::BNGConfig bng_config;
  BNG::BNGEngine bng_engine(bng_config);

  // parse cplx string
  string bngl_str = to_bngl_str_w_custom_orientation(true);
  BNG::Cplx cplx_inst(&bng_engine.get_data());
  int num_errors = BNG::parse_single_cplx_string(
      bngl_str, bng_engine.get_data(),
      cplx_inst
  );
  if (num_errors != 0) {
    throw RuntimeError("While creating canonical name for a Complex, could not parse '" + bngl_str + "'.");
  }
  assert(!cplx_inst.elem_mols.empty());

  // create canonical name
  // we don't care that the complex may not be fully qualified, we cannot know this at this point
  BNG::Species new_species(cplx_inst, bng_engine.get_data(), bng_engine.get_config(), false);
  // order components by name during canonicalization
  new_species.canonicalize(true);

  canonical_name = new_species.name;

  // conversion to Species caused compartment and orientation to be ignored
  if (orientation == Orientation::UP) {
    canonical_name += "'";
  }
  else if (orientation == Orientation::DOWN) {
    canonical_name += ",";
  }

  // compartment_name may not be set even when elementary molecules have their compartments
  if (elementary_molecules.size() != 1 && is_set(compartment_name) && compartment_name != get_primary_compartment_name()) {
    canonical_name = "@" + compartment_name + ":" + canonical_name;
  }

  cached_data_are_uptodate = true;
  return canonical_name;
}


bool Complex::__eq__(const Complex& other) const {

  // cannot use eq_nonarray_attributes here because we don't care about name
  if (orientation != other.orientation) {
    return false;
  }

  return get_canonical_name() == other.get_canonical_name();
}


std::string Complex::to_bngl_str_w_custom_orientation(
    const bool replace_orientation_w_up_down_compartments, const bool ignore_orientation_and_compartment) const {
  string res;
  bool orientation_replaced = false;

  // individual compartments are printed if all elem mols do not use the same compartment
  set<string> used_compartments;
  size_t num_specific_compartments = 0;
  for (const auto& em: elementary_molecules) {
    if (is_set(em->compartment_name)) {
      num_specific_compartments++;
      used_compartments.insert(em->compartment_name);
    }
  }
  bool print_individual_compartments =
      used_compartments.size() > 1 || num_specific_compartments != elementary_molecules.size();

  for (size_t i = 0; i < elementary_molecules.size(); i++) {
    res += elementary_molecules[i]->to_bngl_str(print_individual_compartments);
    if (i + 1 != elementary_molecules.size()) {
      res += ".";
    }
  }

  if (!ignore_orientation_and_compartment) {
    if (!replace_orientation_w_up_down_compartments) {
      if (orientation == Orientation::UP) {
        res += "'";
      }
      else if (orientation == Orientation::DOWN) {
        res += ",";
      }
    }
    else {
      if (orientation == Orientation::UP) {
        res += BNG::MCELL_COMPARTMENT_UP;
        orientation_replaced = true;
      }
      else if (orientation == Orientation::DOWN) {
        res += BNG::MCELL_COMPARTMENT_DOWN;
        orientation_replaced = true;
      }
    }
  }


  if (!print_individual_compartments && used_compartments.size() == 1) {
    string at_str = (orientation_replaced) ? "" : "@";
    if (elementary_molecules.size() != 1) {
      res = at_str + *used_compartments.begin() + ":" + res;
    }
    else {
      res += at_str + *used_compartments.begin();
    }
  }

  return res;
}


bool Complex::is_surf() const {
  for (auto em: elementary_molecules) {
    if (is_set(em->elementary_molecule_type->diffusion_constant_2d)) {
      return true;
    }
  }
  return false;
}


void Complex::check_and_set_compartments() {
  // called after all initialization including parsing of the BNGL name
  // does not check that compartment types are correctly used (surf/vol) because we do not know them yet
  // TODOCOMP: check this in conversion
  if (is_set(compartment_name)) {
    // set the name to all contained elementary molecules
    for (auto& em: elementary_molecules) {
      if (is_set(em->compartment_name) && em->compartment_name != compartment_name) {
        throw ValueError("Cannot override compartment of elementary molecule type of " + em->elementary_molecule_type->name +
            " that uses compartment " + em->compartment_name + " with compartment " + compartment_name +
            " of the Complex " + name + ". If different compartments for individual elementary molecules are needed, " +
            " do not set the Complex's " + NAME_COMPARTMENT_NAME + ".");
      }
      else {
        em->compartment_name = compartment_name;
      }
    }
  }
}


bool Complex::is_species_object() const {
  return dynamic_cast<const Species*>(this) != nullptr;
}


std::shared_ptr<Species> Complex::as_species() {
  return std::make_shared<Species>(*this);
}


const std::string& Complex::get_primary_compartment_name() const {
  if (is_set(compartment_name)) {
    return compartment_name;
  }
  if (elementary_molecules.size() == 1) {
    return elementary_molecules[0]->compartment_name;
  }
  else if (is_surf()) {
    // get the first surface elem mol
    for (const auto& em: elementary_molecules) {
      if (em->is_surf()) {
        return em->compartment_name;
      }
    }
    // unreachable
    throw RuntimeError("Internal error: surface complex " + to_bngl_str() + " does not contain any surface elementary molecules.");
  }
  else {
    // all are volume and must have the same compartment
    return elementary_molecules[0]->compartment_name;
  }
}

} // namespace API
} // namespace MCell
