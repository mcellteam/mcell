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


std::shared_ptr<API::Complex> Complex::construct_from_bng_cplx(
    const BNG::BNGData& bng_data,
    const BNG::Cplx& bng_cplx) {

  std::shared_ptr<API::Complex> res_cplx_inst = API::Complex::construct_empty();

  // convert each molecule instance
  for (const BNG::ElemMol& bmg_mi: bng_cplx.elem_mols) {

    // find molecule type and create an instance
    const BNG::ElemMolType& emt = bng_data.get_elem_mol_type(bmg_mi.elem_mol_type_id);
    std::shared_ptr<ElementaryMoleculeType> api_emt =
        ElementaryMoleculeType::construct_from_bng_elem_mol_type(bng_data, emt);
    assert(is_set(api_emt));

    // prepare a vector of component instances with their bonds set
    std::vector<std::shared_ptr<API::Component>> api_comp_instances;
    for (const BNG::Component& bng_ci: bmg_mi.components) {
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

    // and append instantiated elementary molecule type
    res_cplx_inst->elementary_molecules.push_back(api_emt->inst(api_comp_instances));
  }

  // set compartment
  if (bng_cplx.has_compartment()) {
    if (bng_cplx.has_compartment_class_in_out()) {
      res_cplx_inst->compartment_name = BNG::compartment_id_to_str(bng_cplx.get_compartment_id());
    }
    else {
      res_cplx_inst->compartment_name = bng_data.get_compartment(bng_cplx.get_compartment_id()).name;
    }
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

  // set compartment_name and check that there is only one
  if (is_set(name)) {
    std::vector<std::string> compartments;
    // TODO: we do not want to use the BNG parser at this point or do we?
    get_compartment_names(name, compartments);
    if (!compartments.empty()) {
      if (is_set(compartment_name)) {
        throw ValueError("Complex " + name + " is defined with both compartment name in " +
            NAME_NAME + " and in " + NAME_COMPARTMENT_NAME + ", only one is allowed.");
      }

      std::string single_compartment_name = compartments[0];
      for (size_t i = 1; i < compartments.size(); i++) {
        if (single_compartment_name != compartments[i]) {
          throw ValueError("Complex cannot be in multiple compartments, error for " + name + ".");
        }
      }
      compartment_name = single_compartment_name;
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

    // and copy resulting elementary molecules array
    elementary_molecules = converted_complex->elementary_molecules;
    assert(!is_set(converted_complex->compartment_name) || compartment_name == converted_complex->compartment_name);
  }

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
  if (is_set(compartment_name)) {
    if (name.find('.') == string::npos) {
      canonical_name += "@" + compartment_name;
    }
    else {
      canonical_name = "@" + compartment_name + ":" + canonical_name;
    }
  }

  cached_data_are_uptodate = true;
  return canonical_name;
}


bool Complex::__eq__(const Complex& other) const {

  // cannot use eq_nonarray_attributes here because we don't care about name
  if (orientation != other.orientation ||
      compartment_name != other.compartment_name) {
    return false;
  }

  return get_canonical_name() == other.get_canonical_name();
}


std::string Complex::to_bngl_str_w_custom_orientation(
    const bool replace_orientation_w_up_down_compartments, const bool ignore_orientation_and_compartment) const {
  string res;
  bool add_compartment = false;
  bool orientation_replaced = false;
  if (is_set(name)) {
    res = name;
    if (!ignore_orientation_and_compartment && is_set(compartment_name) && name.find('@') == string::npos) {
      add_compartment = true;
    }
  }
  else {
    for (size_t i = 0; i < elementary_molecules.size(); i++) {
      res += elementary_molecules[i]->to_bngl_str();
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

      if (is_set(compartment_name)) {
        add_compartment = true;
      }
    }
  }

  if (add_compartment) {
    string at_str = (orientation_replaced) ? "" : "@";
    if (name.find('.') == string::npos) {
      res += at_str + compartment_name;
    }
    else {
      res = at_str + compartment_name + ":" + res;
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


bool Complex::is_species_object() const {
  return dynamic_cast<const Species*>(this) != nullptr;
}


std::shared_ptr<Species> Complex::as_species() {
  return std::make_shared<Species>(*this);
}


} // namespace API
} // namespace MCell
