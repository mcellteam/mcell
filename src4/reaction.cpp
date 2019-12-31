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

#include <iostream>

#include "reaction.h"

#include "species_info.h"

using namespace std;

namespace MCell {

void SpeciesWithOrientation::dump_array(const std::vector<SpeciesWithOrientation>& vec, const string ind) {
  for (size_t i = 0; i < vec.size(); i++) {
    cout << ind << i << ": species_id: " << vec[i].species_id << ", orientation:" << vec[i].orientation << "\n";
  }
}

// create mapping for cases when one of the reactants is unchanged in the reaction
void Reaction::update_equivalent_product_indices() {
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
}


uint Reaction::get_num_surf_products(const SpeciesInfo& all_species) const {
  uint res = 0;
  for (const SpeciesWithOrientation& prod: products) {
    if (all_species.get(prod.species_id).is_surf()) {
      res++;
    }
  }
  return res;
}

void Reaction::dump_array(const vector<Reaction>& vec) {
  cout << "Reaction array: " << (vec.empty() ? "EMPTY" : "") << "\n";

  for (size_t i = 0; i < vec.size(); i++) {
    cout << i << ":\n";
    vec[i].dump("  ");
  }
}

void Reaction::dump(const string ind) const {
  cout << ind << "name: \t\t" << name << " [string] \t\t\n";
  cout << ind << "rate_constant: \t\t" << rate_constant << " [float_t] \t\t\n";
  cout << ind << "max_fixed_p: \t\t" << max_fixed_p << " [float_t] \t\t\n";
  cout << ind << "min_noreaction_p: \t\t" << min_noreaction_p << " [float_t] \t\t\n";

  cout << ind << "rectants:\n";
  SpeciesWithOrientation::dump_array(reactants, ind + "  ");

  cout << ind << "products:\n";
  SpeciesWithOrientation::dump_array(products, ind + "  ");
}


} // namespace mcell
