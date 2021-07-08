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

#include <api/component.h>
#include "subsystem.h"

#include "bng/bng.h"

#include "api/component_type.h"
#include "api/elementary_molecule_type.h"
#include "api/elementary_molecule.h"
#include "api/complex.h"

#include "api/api_utils.h"

using namespace std;

namespace MCell {
namespace API {

void Subsystem::dump() const {
  std::cout << to_str() << "\n";
}

std::shared_ptr<Species> Subsystem::get_species_with_id(const species_id_t id) {
  // not very efficient, we may need some caching/map later
  for (auto& s: species) {
    if (s->species_id == id) {
      return s;
    }
  }
  return std::shared_ptr<Species>(nullptr);
}


void Subsystem::unify_and_register_elementary_molecule_types() {
  // then go through Species
  for (std::shared_ptr<Species> s: species) {
    for (size_t i = 0; i < s->elementary_molecules.size(); i++) {
      // note the reference, we might be changing this spared ptr
      std::shared_ptr<ElementaryMoleculeType>& emt = s->elementary_molecules[i]->elementary_molecule_type;

      // do we have this exact object already?
      auto it_ptr_eq = std::find_if(
          elementary_molecule_types.begin(), elementary_molecule_types.end(),
          [&emt] (const std::shared_ptr<ElementaryMoleculeType> emt_existing) {
            return emt_existing.get() == emt.get(); // comparing pointers
          }
      );
      if (it_ptr_eq != elementary_molecule_types.end()) {
        // same object exists
        continue;
      }

      // do we have an object with the same contents already?
      auto it_contents_eq = std::find_if(
          elementary_molecule_types.begin(), elementary_molecule_types.end(),
          [&emt] (const std::shared_ptr<ElementaryMoleculeType> emt_existing) {
            return *emt_existing == *emt; // comparing contents
          }
      );
      if (it_contents_eq != elementary_molecule_types.end()) {
        // same object exists, we must use this one and let shared_ptr ref counting
        // destroy the previous one (creation of the EMT object to be destroyed
        // occurs mainly when defining simple species, so the link from species to this emt
        // is internal and not visible to the user)
        emt = *it_contents_eq;
        continue;
      }

      // one more option is that some of the elementary molecule types that we encounter are not
      // fully initialized, for now expecting that the the fully initialized one is already present
      auto it_name_eq = std::find_if(
          elementary_molecule_types.begin(), elementary_molecule_types.end(),
          [&emt] (const std::shared_ptr<ElementaryMoleculeType> emt_existing) {
            return emt_existing->name == emt->name; // comparing contents
          }
      );
      if (it_name_eq != elementary_molecule_types.end()) {
        // use the one that is initialized
        if (!(*it_name_eq)->all_numerical_attributes_are_unset() &&
            emt->all_numerical_attributes_are_unset()
        ) {
          // update the pointer used in the Species object
          emt = *it_name_eq;
          continue;
        }

        if ((*it_name_eq)->all_numerical_attributes_are_unset() &&
            !emt->all_numerical_attributes_are_unset()
        ) {
          // update the one in elementary_molecule_types using the Species object
          *it_name_eq = emt;
          continue;
        }

        // other cases are error and will be reported in add_elementary_molecule_type
      }

      // no such object is in the emt list, add it (reports error if such object already exists)
      add_elementary_molecule_type(emt);
    }
  }
}


void Subsystem::load_bngl_molecule_types_and_reaction_rules(
    const std::string& file_name,
    const std::map<std::string, double>& parameter_overrides) {

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
  for (const BNG::ElemMolType& mt: bng_data.get_elem_mol_types()) {
    auto res_mt = ElementaryMoleculeType::construct_from_bng_elem_mol_type(bng_data, mt);
    append_to_vec_canonical_name(elementary_molecule_types, res_mt);
  }

  // reaction rules
  for (const BNG::RxnRule& rr: bng_data.get_rxn_rules()) {
    convert_reaction_rule(bng_data, rr);
  }
}


void Subsystem::convert_reaction_rule(const BNG::BNGData& bng_data, const BNG::RxnRule& bng_rr) {

  // always a standard rxn
  if (bng_rr.type != BNG::RxnType::Standard) {
    throw RuntimeError("Unexpected type of reaction from BNGL file for reaction " +  bng_rr.name + ".");
  }

  auto res_rr = make_shared<API::ReactionRule>(DefaultCtorArgType());
  if (bng_rr.name != "") {
    res_rr->name = bng_rr.name;
  }

  // RxnRule
  res_rr->fwd_rate = bng_rr.base_rate_constant;

  for (const BNG::Cplx& inst: bng_rr.reactants) {
    // MCell3R accepts reactants with any orientation
    res_rr->reactants.push_back(
        Complex::construct_from_bng_cplx_w_orientation(bng_data, inst, Orientation::ANY));
  }

  for (const BNG::Cplx& inst: bng_rr.products) {
    // MCell3R always creates products with the orientation up
    res_rr->products.push_back(
        Complex::construct_from_bng_cplx_w_orientation(bng_data, inst, Orientation::UP));
  }

  // allow reactions with identical names
  append_to_vec_canonical_name(reaction_rules, res_rr);
}


} // namespace API
} // namespace MCell

