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
#include "api/elementary_molecule.h"
#include "api/elementary_molecule_type.h"
#include "bng/bng.h"

using namespace std;

namespace MCell {
namespace API {

const std::string& Complex::get_canonical_name() const {
  if (cached_data_are_uptodate && canonical_name != "") {
    return canonical_name;
  }

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
    const bool replace_orientation_w_up_down_compartments, const bool ignore_orientation) const {
  string res;
  bool add_compartment = false;
  bool orientation_replaced = false;
  if (is_set(name)) {
    res = name;
    if (is_set(compartment_name) && name.find('@') == string::npos) {
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
