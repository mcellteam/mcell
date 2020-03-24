/******************************************************************************
 *
 * Copyright (C) 2019 by
 * The Salk Institute for Biological Studies and
 * Pittsburgh Supercomputing Center, Carnegie Mellon University
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

#include "reaction.h"

#include <iostream>

#include "species.h"
#include "partition.h"

using namespace std;

namespace MCell {

void SpeciesWithOrientation::dump_array(const std::vector<SpeciesWithOrientation>& vec, const string ind) {
  for (size_t i = 0; i < vec.size(); i++) {
    cout << ind << i << ": species_id: " << vec[i].species_id << ", orientation:" << vec[i].orientation << "\n";
  }
}


void Rxn::initialize(const RxnClass& reaction) {
  debug_check_reactants_against_reaction(reaction);

  update_equivalent_product_indices();

  // MCELL3 compatibility - reorder products so that case such as
  // CaM -> Ca + CaM becomes CaM -> CaM + Ca
  move_reused_reactants_to_be_the_first_products();
}

// asserts in debug mode if the reactants are different
void Rxn::debug_check_reactants_against_reaction(const RxnClass& reaction) {
  assert(false && "BNGTODO");
  /*
  assert(reaction.reactants.size() == reactants.size());
  assert(reactants.size() >= 1 && reactants.size() <= 2);
  if (reactants.size() == 1) {
    assert(reactants[0] == reaction.reactants[0]);
  }
  */
}

uint Rxn::get_num_surf_products(const SpeciesInfo& all_species) const {
  assert(false && "BNGTODO");
  /*
  uint res = 0;
  for (const SpeciesWithOrientation& prod: products) {
    if (all_species.get(prod.species_id).is_surf()) {
      res++;
    }
  }
  return res;*/
}

// create mapping for cases when one of the reactants is unchanged in the reaction
void Rxn::update_equivalent_product_indices() {
  assert(false && "BNGTODO");
  /*
  for (SpeciesWithOrientation& product: products) {
    product.equivalent_product_or_reactant_index = INDEX_INVALID;
  }

  for (uint ri = 0; ri < reactants.size(); ri++) {
    reactants[ri].equivalent_product_or_reactant_index = INDEX_INVALID;

    for (uint pi = 0; pi < products.size(); pi++) {
      if (reactants[ri].equivalent_product_or_reactant_index == INDEX_INVALID &&
          products[pi].equivalent_product_or_reactant_index == INDEX_INVALID &&
          reactants[ri] == products[pi]) {

        reactants[ri].equivalent_product_or_reactant_index = pi;
        products[pi].equivalent_product_or_reactant_index = ri;
      }
    }
  }
  */
}


void Rxn::move_reused_reactants_to_be_the_first_products() {
  assert(false && "BNGTODO");
  /*
  // for each reactant (from the end since we want the products to be ordered in the same way)
  for (int ri = reactants.size() - 1; ri >= 0; ri--) {
    if (reactants[ri].equivalent_product_or_reactant_index != INDEX_INVALID) {
      // move product to the front
      uint index = reactants[ri].equivalent_product_or_reactant_index;
      SpeciesWithOrientation prod = products[index];
      products.erase(products.begin() + index);
      products.insert(products.begin(), prod);

      // update mapping (inefficient, but used only in initialization)
      update_equivalent_product_indices();
    }
  }
  */
}

void Rxn::dump(const string ind) const {
  assert(false && "TODO");
  /*
  cout << ind << "name: \t\t" << name << " [string] \t\t\n";
  cout << ind << "rate_constant: \t\t" << rate_constant << " [float_t] \t\t\n";

  cout << ind << "rectants:\n";
  SpeciesWithOrientation::dump_array(reactants, ind + "  ");

  cout << ind << "products:\n";
  SpeciesWithOrientation::dump_array(products, ind + "  ");
  */
}


static std::string pathway_players_to_string(
    const Partition& p, const std::vector<SpeciesWithOrientation>& players
) {
  string res;

  for (uint i = 0; i < players.size(); i++) {
    res += p.all_species.get(players[i].species_id).name;
    if (i != players.size() - 1) {
      res += " + ";
    }
  }
  return res;
}


std::string Rxn::to_string(const Partition& p) const {
  assert(false && "TODO");
  //return pathway_players_to_string(p, reactants) + " -> " + pathway_players_to_string(p, products);
}

void RxnClass::dump_array(const vector<RxnClass>& vec) {
  cout << "Reaction array: " << (vec.empty() ? "EMPTY" : "") << "\n";

  for (size_t i = 0; i < vec.size(); i++) {
    cout << i << ":\n";
    vec[i].dump("  ");
  }
}

void RxnClass::dump(const string ind) const {
  cout << ind << "max_fixed_p: \t\t" << max_fixed_p << " [float_t] \t\t\n";
  cout << ind << "min_noreaction_p: \t\t" << min_noreaction_p << " [float_t] \t\t\n";
  SpeciesWithOrientation::dump_array(reactants, ind + "  ");

  for (const Rxn& rxn: reactions) {
    rxn.dump(ind);
  }
}




} // namespace mcell
