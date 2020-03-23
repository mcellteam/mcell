/*
 * rule.cpp
 *
 *  Created on: Jan 9, 2020
 *      Author: ahusar
 */

#include <iostream>
#include <sstream>

#include "rxn_rule.h"

using namespace std;

namespace BNG {



CplxMolIndex* RxnRule::get_assigned_reactant_for_product(CplxMolIndex& product_cmi) {

}


// score:
//  - start with -1
//  - components are the same: +1
//  - each fitting state: +1
void RxnRule::find_most_fitting_unassigned_product(CplxMolIndex& reactant_cmi) {
  const MolInstance& mi = get_reactant_mol(reactant_cmi);

  int best_score = -1;
  CplxMolIndex best_cmi;

  for (uint complex_index = 0; complex_index < products.size(); complex_index++) {
    for (uint molecule_index = 0; molecule_index < products[complex_index].mol_patterns.size(); molecule_index++) {

      CplxMolIndex product_cmi(complex_index, molecule_index);

      // skip if already assigned
      CplxMolIndex* found_cmi = get_assigned_reactant_for_product(product_cmi);
      if (found_cmi != nullptr) {
        continue;
      }


      int score = -1;

      //mapping

    }
  }
}



// check if it makes sense to compute mapping at all
bool RxnRule::has_same_molecules_in_reactants_and_products() {
  map<mol_type_id_t, int> reactant_types, product_types;

  for (const CplxInstance& ci: reactants) {
    for (const MolInstance& mi: ci.mol_patterns) {
      if (reactant_types.count(mi.mol_type_id)) {
        reactant_types[mi.mol_type_id]++;
      }
      else {
        reactant_types[mi.mol_type_id] = 0;
      }
    }
  }

  for (const CplxInstance& ci: products) {
    for (const MolInstance& mi: ci.mol_patterns) {
      if (reactant_types.count(mi.mol_type_id)) {
        product_types[mi.mol_type_id]++;
      }
      else {
        product_types[mi.mol_type_id] = 0;
      }
    }
  }

  return reactant_types == product_types;
}



bool RxnRule::compute_reactants_products_mapping(std::stringstream& msgs) {

  if (!has_same_molecules_in_reactants_and_products()) {
    mol_instances_are_maintained = false;
    return false;
  }


  for (uint complex_index = 0; complex_index < reactants.size(); complex_index++) {
    for (uint molecule_index = 0; molecule_index < reactants[complex_index].mol_patterns.size(); molecule_index++) {

      find_most_fitting_unassigned_product(CplxMolIndex(complex_index, molecule_index));
    }
  }

  return true;
}


void RxnRule::dump_complex_instance_vector(const BNGData& bng_data, const ComplexInstanceVector& complexes) const {

  for (size_t i = 0; i < complexes.size(); i++) {
    complexes[i].dump(bng_data);

    if (i != complexes.size() - 1) {
      cout << " + ";
    }
  }
}


void RxnRule::dump(const BNGData& bng_data) const {
  if (name != "") {
    cout << name << " ";
  }
  dump_complex_instance_vector(bng_data, reactants);

  cout << " -> ";
  dump_complex_instance_vector(bng_data, products);

  cout << " " << rxn_rate;
}

} /* namespace BNG */
