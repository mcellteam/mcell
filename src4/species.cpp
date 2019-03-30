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

#include "species.h"

using namespace std;

namespace mcell {

void species_t::dump(const string ind) const {
  cout << ind <<"species_id: \t\t" << species_id << " [uint16_t] \t\t/* Unique ID for this species */\n";
  cout << ind <<"mcell_species_id: \t\t" << mcell3_species_id << " [uint32_t] \t\t/* Unique ID for this species from mcell3 representation*/\n";
  cout << ind <<"name: *\t\t" << name << " [string] \t\t/* Symbol table entry (name) */\n";
  cout << ind <<"D: \t\t" << D << " [float_t] \t\t/* Diffusion constant */\n";
  cout << ind <<"space_step: \t\t" << space_step << " [float_t] \t\t/* Characteristic step length */\n";
  cout << ind <<"time_step: \t\t" << time_step << " [float_t] \t\t/* Minimum (maximum?) sensible timestep */\n";
}

void species_t::dump_array(const vector<species_t>& vec) {
	cout << "Species array: " << (vec.empty() ? "EMPTY" : "") << "\n";

	for (size_t i = 0; i < vec.size(); i++) {
		cout << i << ":\n";
		vec[i].dump("  ");
	}
}

} // namespace mcell
