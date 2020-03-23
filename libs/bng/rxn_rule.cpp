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



bool RxnRule::find_assigned_reactant_for_product(const CplxMolIndex& product_cmi, CplxMolIndex& reactant_cmi) const {
  // this is not a time critical search
  for (const CMIndexPair& cmi_pair: mapping) {
    if (product_cmi == cmi_pair.product_cmi) {
      reactant_cmi = cmi_pair.reactant_cmi;
      return true;
    }
  }
  return false;
}


// Finds a matching already not assigned product,
// the product must be
// NOTE: BNGL2.pl provides more detailed reporting, see tests N220, N230-N232
bool RxnRule::find_most_fitting_unassigned_product(const CplxMolIndex& reactant_cmi, CplxMolIndex& best_product_cmi) const {
  const MolInstance& reactant_mol_inst = get_reactant_mol_inst(reactant_cmi);

  int best_score = -1;
  CplxMolIndex best_cmi;

  for (uint complex_index = 0; complex_index < products.size(); complex_index++) {
    for (uint molecule_index = 0; molecule_index < products[complex_index].mol_patterns.size(); molecule_index++) {

      CplxMolIndex product_cmi(complex_index, molecule_index);

      // must have the same molecule type
      const MolInstance& product_mol_inst = get_product_mol_inst(product_cmi);
      if (reactant_mol_inst.mol_type_id != product_mol_inst.mol_type_id) {
        continue;
      }

      // and must not be assigned
      CplxMolIndex found_cmi_ignored;
      bool found = find_assigned_reactant_for_product(product_cmi, found_cmi_ignored);
      if (found) {
        continue;
      }

      // ok, this product was not mapped yet
      int num_explicitly_listed_components_in_reactant = 0;
      int num_same_explicitly_listed_components = 0;
      int num_same_component_states = 0; // when the components are expl. listed

      assert(reactant_mol_inst.component_instances.size() == product_mol_inst.component_instances.size());
      for (uint i = 0; i < reactant_mol_inst.component_instances.size(); i++) {

        const ComponentInstance& reactant_comp_inst = reactant_mol_inst.component_instances[i];
        const ComponentInstance& product_comp_inst = product_mol_inst.component_instances[i];

        // component must be explicitly listed on both sides to be considered
        if (reactant_comp_inst.explicitly_listed_in_pattern) {
          num_explicitly_listed_components_in_reactant++;

          // the same component must be explicitly listed and
          // if state is specified, it must be set on both sides
          if (product_comp_inst.explicitly_listed_in_pattern &&
              reactant_comp_inst.state_is_set() == product_comp_inst.state_is_set())  {

            num_same_explicitly_listed_components++;
            if (reactant_comp_inst.state_id == product_comp_inst.state_id) {
              num_same_component_states++;
            }
          }
        }
        else if (product_comp_inst.explicitly_listed_in_pattern) {
          // component is not listed in reactants, just in products
          num_same_explicitly_listed_components--;
        }
      }

      // all listed components match?
      if (num_explicitly_listed_components_in_reactant == num_same_explicitly_listed_components) {
        if (num_same_component_states > best_score) {
          best_score = num_same_component_states;
          best_cmi = product_cmi;
        }
      }
    }
  }

  if (best_score != -1) {
    best_product_cmi = best_cmi;
    return true;
  }
  else {
    return false;
  }
}



// check if it makes sense to compute mapping at all
bool RxnRule::has_same_molecules_in_reactants_and_products() const {
  map<mol_type_id_t, int> reactant_types, product_types;

  for (const CplxInstance& ci: reactants) {
    for (const MolInstance& mi: ci.mol_patterns) {
      if (reactant_types.count(mi.mol_type_id)) {
        reactant_types[mi.mol_type_id]++;
      }
      else {
        reactant_types[mi.mol_type_id] = 1;
      }
    }
  }

  for (const CplxInstance& ci: products) {
    for (const MolInstance& mi: ci.mol_patterns) {
      if (product_types.count(mi.mol_type_id)) {
        product_types[mi.mol_type_id]++;
      }
      else {
        product_types[mi.mol_type_id] = 1;
      }
    }
  }

  return reactant_types == product_types;
}



bool RxnRule::compute_reactants_products_mapping(const BNGData& bng_data, std::ostream& out) {

  if (!has_same_molecules_in_reactants_and_products()) {
    mol_instances_are_maintained = false;
    return true; // this is ok because molecules are transformed in this rule
  }


  for (uint complex_index = 0; complex_index < reactants.size(); complex_index++) {
    for (uint molecule_index = 0; molecule_index < reactants[complex_index].mol_patterns.size(); molecule_index++) {

      CplxMolIndex reactant_cmi = CplxMolIndex(complex_index, molecule_index);
      CplxMolIndex product_cmi;
      bool found = find_most_fitting_unassigned_product(reactant_cmi, product_cmi);

      if (found) {
        mapping.push_back(CMIndexPair(reactant_cmi, product_cmi));
      }
      else {
        const MolInstance& mol_inst = get_reactant_mol_inst(reactant_cmi);

        out << "Did not find a matching molecule in products for reactant molecule ";
        mol_inst.dump(bng_data, true, out);
        out << " listed as complex " << reactant_cmi.cplx_index << " and molecule " << reactant_cmi.mol_index << ".";

        return false;
      }

    }
  }

  return true;
}


void RxnRule::dump_complex_instance_vector(const BNGData& bng_data, const CplxInstanceVector& complexes) const {

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
