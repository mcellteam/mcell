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

#include "molecule.h"

using namespace std;

namespace mcell {


void base_molecule_t::dump_base(const std::string ind) const {
  cout << ind <<"flags: \t\t" << flags << "[uint16_t]\n";
  cout << ind <<"species_id: \t\t" << species_id << " [species_id_t]\n";
}

void volume_molecule_t::dump(const std::string ind) const {
  cout << ind <<"pos: \t\t" << pos << "[vec3_t]\n";
  cout << ind <<"subpartition_index: \t\t" << subpartition_index << " [uint32_t]\n";
  dump_base(ind);
}


} /* namespace mcell */
