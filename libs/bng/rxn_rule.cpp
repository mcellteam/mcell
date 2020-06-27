/*
 * rule.cpp
 *
 *  Created on: Jan 9, 2020
 *      Author: ahusar
 */

#include <iostream>
#include <sstream>

#include "bng/rxn_rule.h"
#include "bng/rxn_class.h"

#include "bng/species.h"
#include "bng/species_container.h"

#include "bng/graph.h"

using namespace std;

namespace BNG {

void RxnRule::finalize() {
  assert(id != RXN_RULE_ID_INVALID);

  bool simple = true;

  // finalize all reactants and products
  for (CplxInstance& ci: reactants) {
    ci.finalize();
    simple &= ci.is_simple();
  }

  num_surf_products = 0;
  for (CplxInstance& ci: products) {
    ci.finalize();
    simple &= ci.is_simple();
    if (ci.is_surf()) {
      num_surf_products++;
    }
  }

  if (simple) {
    set_flag(RXN_FLAG_SIMPLE);
  }

  compute_reactants_products_mapping();

  // for MCell3 compatibility
  move_products_that_are_also_reactants_to_be_the_first_products();

  set_finalized();
}


static void merge_graphs(Graph& srcdst, const Graph& src) {
  assert(false && "TODO");
}


void RxnRule::create_products_for_complex_rxn(
    const std::vector<const CplxInstance*>& input_reactants,
    std::vector<CplxInstance>& created_products
) const {

  assert(mol_instances_are_fully_maintained && "Assuming this for now");

  assert(input_reactants.size() == reactants.size());
  assert(input_reactants.size() == 1 || input_reactants.size() == 2);

  // merge input reactant graphs
  Graph reactants_graph = input_reactants[0]->get_graph();
  if (input_reactants.size() == 2) {
    merge_graphs(reactants_graph, input_reactants[1]->get_graph());
  }

  // merge reactant reactant patterns
  // this can be precomputed
  Graph pattern_graph = reactants[0].get_graph();
  if (reactants.size() == 2) {
    merge_graphs(pattern_graph, reactants[1].get_graph());
  }

  // compute mapping reactant pattern -> reactant
  MultipleMappingsVector pattern_to_reactants_mappings;
  get_subgraph_isomorphism_mappings(
      pattern_graph, // pattern
      reactants_graph, // actual reactant
      false, // do not stop with first match
      pattern_to_reactants_mappings
  );
  assert(pattern_to_reactants_mappings.size() == 1 &&
      "We do not support multiple matches yet and then there must be at least one match.");


  // create graph from products


  // compute mapping reactant patterns -> product patterns
  // this can be precomputed
  // nodes in the reaction products graph have an extra option that they can ignore the state
  // find all mappings and compute score for each of those mappings so that we get the best match
  // -> each matching state is a positive point

  //find_best_reactants_to_product_mapping()


  // manipulate nodes using information about products
  //
  // for each component in product graph:
  //   find corresponding component in reactant pattern
  //     if there is no such mapping, ignore it because it might be a component of a
  //     new molecule instance
  //   find component corresponding to reactant pattern in reactants_graph
  //     this mapping must exist because we matched the pattern graph to reactants graph
  //   manipulate state and or bond
  //     removal of a bond should be ok
  //
  // for each molecule in product graph:
  //   if there is no mapping to reactant pattern graph:
  //     check that all components that have state are connected and defined in the product graph
  //     create a new molecule instance in the reactants_graph
  //     for each component of the molecule in product graph:
  //       add it and connect to molecule
  //       if it has a bond:
  //         if does the second target already exist:
  //           create a bond
  //


  // and create products, each disconnected graph in the result is a
  // separate complex instance

}


bool RxnRule::is_cplx_reactant_on_both_sides_of_rxn(const uint index) const {
  assert(is_finalized());
  for (const CplxIndexPair& cplx_index_pair: cplx_mapping) {
    if (index == cplx_index_pair.reactant_index) {
      return true;
    }
  }
  return false;
}


bool RxnRule::is_cplx_product_on_both_sides_of_rxn(const uint index) const {
  assert(is_finalized());
  for (const CplxIndexPair& cplx_index_pair: cplx_mapping) {
    if (index == cplx_index_pair.product_index) {
      return true;
    }
  }
  return false;
}


bool RxnRule::find_assigned_mol_reactant_for_product(const CplxMolIndex& product_cmi, CplxMolIndex& reactant_cmi) const {
  // this is not a time critical search
  for (const CMIndexPair& cmi_pair: mol_mapping) {
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
bool RxnRule::find_most_fitting_unassigned_mol_product(const CplxMolIndex& reactant_cmi, CplxMolIndex& best_product_cmi) const {
  const MolInstance& reactant_mol_inst = get_mol_reactant(reactant_cmi);

  int best_score = -1;
  CplxMolIndex best_cmi;

  for (uint complex_index = 0; complex_index < products.size(); complex_index++) {
    for (uint molecule_index = 0; molecule_index < products[complex_index].mol_instances.size(); molecule_index++) {

      CplxMolIndex product_cmi(complex_index, molecule_index);

      // must have the same molecule type
      const MolInstance& product_mol_inst = get_mol_product(product_cmi);
      if (reactant_mol_inst.mol_type_id != product_mol_inst.mol_type_id) {
        continue;
      }

      // and must not be assigned
      CplxMolIndex found_cmi_ignored;
      bool found = find_assigned_mol_reactant_for_product(product_cmi, found_cmi_ignored);
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


// check if it makes sense to compute molecule_mapping at all
bool RxnRule::has_same_mols_in_reactants_and_products() const {
  map<mol_type_id_t, int> reactant_types, product_types;

  for (const CplxInstance& ci: reactants) {
    for (const MolInstance& mi: ci.mol_instances) {
      if (reactant_types.count(mi.mol_type_id)) {
        reactant_types[mi.mol_type_id]++;
      }
      else {
        reactant_types[mi.mol_type_id] = 1;
      }
    }
  }

  for (const CplxInstance& ci: products) {
    for (const MolInstance& mi: ci.mol_instances) {
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


bool RxnRule::find_assigned_cplx_reactant_for_product(const uint product_index, uint& reactant_index) const {
  // this is not a time critical search
  for (const CplxIndexPair& cplx_index_pair: cplx_mapping) {
    if (product_index == cplx_index_pair.product_index) {
      reactant_index = cplx_index_pair.reactant_index;
      return true;
    }
  }
  return false;
}


// a matching reactant and product must be identical
void RxnRule::compute_cplx_reactants_products_mapping() {

  cplx_mapping.clear();

  for (uint ri = 0; ri < reactants.size(); ri++) {
    for (uint pi = 0; pi < products.size(); pi++) {
      uint index_ignored;
      if (find_assigned_cplx_reactant_for_product(pi, index_ignored)) {
        // already used
        continue;
      }

      if (reactants[ri].matches_fully(products[pi], true)) {
        cplx_mapping.push_back(CplxIndexPair(ri, pi));
        // reactant was mapped, continue with the next reactant
        break;
      }
    }
  }
}


bool RxnRule::compute_mol_reactants_products_mapping(MolInstance& not_matching_mol_inst, CplxMolIndex& not_matching_cmi) {
  mol_mapping.clear();

  mol_instances_are_fully_maintained = has_same_mols_in_reactants_and_products();

  for (uint complex_index = 0; complex_index < reactants.size(); complex_index++) {
    for (uint molecule_index = 0; molecule_index < reactants[complex_index].mol_instances.size(); molecule_index++) {

      CplxMolIndex reactant_cmi = CplxMolIndex(complex_index, molecule_index);
      CplxMolIndex product_cmi;
      bool found = find_most_fitting_unassigned_mol_product(reactant_cmi, product_cmi);

      if (found) {
        mol_mapping.push_back(CMIndexPair(reactant_cmi, product_cmi));
      }
      else if (mol_instances_are_fully_maintained) {
        // reporting error only if there should be a full match
        not_matching_mol_inst = get_mol_reactant(reactant_cmi);
        not_matching_cmi = reactant_cmi;

        return false;
      }
    }
  }

  return true;
}


bool RxnRule::compute_reactants_products_mapping() {

  compute_cplx_reactants_products_mapping();

  MolInstance not_matching_mol_inst_ignored;
  CplxMolIndex not_matching_cmi_ignored;
  bool ok = compute_mol_reactants_products_mapping(not_matching_mol_inst_ignored, not_matching_cmi_ignored);
  assert(ok);
  return ok;
}


bool RxnRule::compute_reactants_products_mapping_w_error_output(const BNGData& bng_data, std::ostream& out) {

  compute_cplx_reactants_products_mapping();

  // NOTE: we might need to direct the molecule mapping using cplx mapping,
  // but let's see later

  MolInstance not_matching_mol_inst;
  CplxMolIndex not_matching_cmi;
  bool ok = compute_mol_reactants_products_mapping(not_matching_mol_inst, not_matching_cmi);
  if (!ok) {
    out << "Did not find a matching molecule in products for reactant molecule ";
    out << not_matching_mol_inst.to_str(bng_data, true);
    out << " listed as complex " << not_matching_cmi.cplx_index << " and molecule " << not_matching_cmi.mol_index << ".";
  }
  return ok;
}


void RxnRule::move_products_that_are_also_reactants_to_be_the_first_products() {

  // for each reactant (from the end since we want the products to be ordered in the same way)
  for (int pi = products.size() - 1; pi > 0; pi--) {
    uint ri;
    bool found = find_assigned_cplx_reactant_for_product(pi, ri);

    if (found) {
      // move product to the front
      CplxInstance prod = products[pi];
      products.erase(products.begin() + pi);
      products.insert(products.begin(), prod);

      // update mapping
      compute_reactants_products_mapping();
    }
  }
}


bool RxnRule::species_can_be_reactant(const species_id_t id, const SpeciesContainer& all_species) {

  // check caches first
  if (species_applicable_as_reactants.count(id) != 0) {
    return true;
  }
  if (species_not_applicable_as_reactants.count(id) != 0) {
    return false;
  }

  // need to find out
  const CplxInstance& inst = all_species.get_as_cplx_instance(id);

  // at least one should match
  bool matches = false;
  for (const CplxInstance& reactant: reactants) {
    if (reactant.matches_pattern(inst, true)) {
      matches = true;
      break;
    }
    else {
      matches = false;
    }
  }

  if (matches) {
    species_applicable_as_reactants.insert_unique(id);
  }
  else {
    species_not_applicable_as_reactants.insert_unique(id);
  }

  return matches;
}


bool RxnRule::species_is_both_bimol_reactants(const species_id_t id, const SpeciesContainer& all_species) {

  if (!is_bimol()) {
    return false;
  }

  // check if the species can be a reactant at all
  if (!species_can_be_reactant(id, all_species)) {
    return false;
  }

  // then the reactants must be identical (this can be precomputed)
  bool res = reactants[0].matches_pattern(reactants[1], true);
  assert(res == reactants[1].matches_pattern(reactants[0], true) && "Pattern identity must be bijective");
  return res;
}


bool RxnRule::update_variable_rxn_rate(const float_t current_time, const RxnClass* requester) {
  if (!may_update_rxn_rate()) {
    return false;
  }
  assert(!variable_rates.empty());
  assert(next_variable_rate_index < variable_rates.size());

  if (variable_rates[next_variable_rate_index].time > current_time
      && !cmp_eq(variable_rates[next_variable_rate_index].time, current_time) ) {
    return false;
  }

  // find which time to use - the highest but still smaller than the following one
  size_t current_index = next_variable_rate_index;
  while (current_index < variable_rates.size() &&
          (current_time > variable_rates[current_index + 1].time ||
           cmp_eq(current_time, variable_rates[current_index + 1].time)
        )
  ) {
    current_index++;
  }

  // should we use the last entry?
  if (current_index == variable_rates.size() &&
      next_variable_rate_index == current_index - 1) {
    current_index = next_variable_rate_index;
  }

  // current_time >= time for next change
  rate_constant = variable_rates[current_index].rate_constant;
  next_variable_rate_index = current_index + 1;

  // notify parents that update is needed
  for (RxnClass* user: rxn_classes_where_used) {
    // do not call update on the class that called us
    if (user != requester) {
      user->update_rxn_rates_if_needed(current_time);
    }
  }

  return true;
}


std::string RxnRule::to_str(const BNGData& bng_data) const {
  stringstream ss;
  ss << name << ": ";

  ss << complex_instance_vector_to_str(bng_data, reactants);
  ss << " -> ";
  ss << complex_instance_vector_to_str(bng_data, products);

  ss << " " << rate_constant;

  return ss.str();
}


std::string RxnRule::complex_instance_vector_to_str(const BNGData& bng_data, const CplxInstanceVector& complexes) const {
  stringstream ss;
  for (size_t i = 0; i < complexes.size(); i++) {
    ss << complexes[i].to_str(bng_data, is_surf_rxn());

    if (i != complexes.size() - 1) {
      ss << " + ";
    }
  }
  return ss.str();
}


std::string RxnRule::reactants_to_str(const BNGData& bng_data) const {
  return complex_instance_vector_to_str(bng_data, reactants);
}


std::string RxnRule::products_to_str(const BNGData& bng_data) const {
  return complex_instance_vector_to_str(bng_data, products);
}


void RxnRule::dump_complex_instance_vector(
    const BNGData& bng_data, const CplxInstanceVector& complexes,
    const std::string ind) const {

  for (size_t i = 0; i < complexes.size(); i++) {
    cout << ind << "CplxInstance " << i << ":\n";
    complexes[i].dump(true, ind + "  ");
    cout << "\n";
  }
}


void RxnRule::dump(const BNGData& bng_data, const bool for_diff, const std::string ind) const {
  if (!for_diff) {
    cout << ind << to_str(bng_data);
  }
  else {
    cout << ind << "name: " << name << "\n";
    cout << ind << "id: " << id << "\n";

    cout << ind << "type: ";
    switch (type) {
      case RxnType::Standard:
        cout << "Standard";
        break;
      case RxnType::AbsorbRegionBorder:
        cout << "AbsorbRegionBorder";
        break;
      case RxnType::Reflect:
        cout << "Reflect";
        break;
      case RxnType::Transparent:
        cout << "Transparent";
        break;
      default:
        assert(false);
    }
    cout << "\n";

    cout << ind << "rate_constant: " << rate_constant << "\n";
    cout << ind << "variable_rates.size: " << variable_rates.size() << "\n";

    if (!variable_rates.empty()) {
      for (size_t i = 0; i < variable_rates.size(); i++) {
        cout << ind + "  " << "t: " << variable_rates[i].time << ", r: " << variable_rates[i].rate_constant << "\n";
      }
    }

    cout << ind << "mol_instances_are_fully_maintained: " << mol_instances_are_fully_maintained << "\n";

    cout << ind << "reactants:\n";
    dump_complex_instance_vector(bng_data, reactants, ind);
    cout << ind << "products:\n";
    dump_complex_instance_vector(bng_data, products, ind);
  }
}


} /* namespace BNG */
