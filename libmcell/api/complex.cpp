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
#include "api/elementary_molecule_instance.h"
#include "api/elementary_molecule_type.h"
#include "bng/bng.h"

using namespace std;

namespace MCell {
namespace API {

void Complex::set_canonical_name_if_needed() const {
  if (cached_data_are_uptodate && canonical_name != "") {
    return;
  }

  // get BNGL string and parse it, then get canonical name
  BNG::BNGConfig bng_config;
  BNG::BNGEngine bng_engine(bng_config);

  // parse cplx string
  BNG::Cplx cplx_inst(&bng_engine.get_data());
  int num_errors = BNG::parse_single_cplx_string(
      to_bngl_str(), bng_engine.get_data(),
      cplx_inst
  );
  if (num_errors != 0) {
    throw RuntimeError("While creating canonical name for a Complex, could not parse '" + to_bngl_str() + "'.");
  }
  assert(!cplx_inst.mol_instances.empty());

  // create canonical name
  // we don't care that the complex may not be fully qualified, we cannot know this at this point
  BNG::Species new_species(cplx_inst, bng_engine.get_data(), bng_engine.get_config(), false);
  // order components by name during canonicalization
  new_species.canonicalize(true);

  canonical_name = new_species.name;
  cached_data_are_uptodate = true;
}


bool Complex::__eq__(const Complex& other) const {

  if (orientation != other.orientation ||
      compartment_name != other.compartment_name) {
    return false;
  }

  set_canonical_name_if_needed();
  other.set_canonical_name_if_needed();

  return canonical_name == other.canonical_name;
}


std::string Complex::to_bngl_str() const {
  if (is_set(name)) {
    return name;
  }
  else {
    std::string res;
    for (size_t i = 0; i < elementary_molecule_instances.size(); i++) {
      res += elementary_molecule_instances[i]->to_bngl_str();
      if (i + 1 != elementary_molecule_instances.size()) {
        res += ".";
      }
    }

    if (orientation == Orientation::UP) {
      res += "'";
    }
    else if (orientation == Orientation::DOWN) {
      res += ",";
    }

    return res;
  }
}


bool Complex::is_surf() const {
  for (auto em: elementary_molecule_instances) {
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
