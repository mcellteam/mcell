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
#include "world_constants.h"

using namespace std;

namespace MCell {

void SpeciesWithOrientation::dump_array(const std::vector<SpeciesWithOrientation>& vec, const string ind) {
  for (size_t i = 0; i < vec.size(); i++) {
    cout << ind << i << ": species_id: " << vec[i].species_id << ", orientation:" << vec[i].orientation << "\n";
  }
}


uint Reaction::get_num_surf_products(const WorldConstants& world_contants) const {
  uint res = 0;
  for (const SpeciesWithOrientation& prod: products) {
    if (world_contants.get_species(prod.species_id).is_surf()) {
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
