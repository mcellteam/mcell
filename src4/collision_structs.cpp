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
#include <string>
#include <sstream>

#include "collision_structs.h"

#include "partition.h"
#include "molecule.h"

using namespace std;

namespace MCell {


void Collision::dump(Partition& p, const std::string ind) const {
  cout << ind << "diffused_molecule:\n";
  p.get_m(diffused_molecule_id).dump(ind + "  ");
  if (type == CollisionType::VOLMOL_VOLMOL) { // dump is rather limited for now, goes not deal with all types
    cout << ind << "colliding_molecule:\n";
    p.get_m(colliding_molecule_id).dump(ind + "  ");
    cout << ind << "reaction:";

    if (rx != nullptr) {
      rx->dump(ind + "  ");
    }
    else if (rxn_class != nullptr) {
      rxn_class->dump(ind + "  ");
    }
  }
  else {
    cout << ind << "colliding_wall_index: " << colliding_wall_index << "\n";
  }

  cout << "time: \t\t" << time << " [float_t] \t\t\n";
  cout << "position: \t\t" << pos << " [vec3_t] \t\t\n";
}

void Collision::dump(
    const Partition& p,
    const std::string extra_comment,
    const uint64_t iteration
) const {

  cout << extra_comment << "it:" << iteration << ", ";

  if (is_mol_mol_reaction()) {
    cout <<
      "bimol rxn" <<
      ", idA:"  << diffused_molecule_id <<
      ", idB:"  << colliding_molecule_id <<
      //TODO in mcell3 ", rxn: " << rx->to_string(p) <<
      ", time: " << time;

      if (type != CollisionType::SURFMOL_SURFMOL) {
        cout << ", pos " << pos;
      }
  }
  else if (is_unimol_reaction()) {
    cout <<
      "unimol rxn" <<
      ", idA:"  << diffused_molecule_id <<
      //TODO in mcell3 ", rxn: " << rx->to_string(p) <<
      ", time: " << time;
  }
  else if (is_wall_collision()) {
    cout << "Unimol (TODO)";
  }
  else {
    assert(false);
  }
  cout << "\n";
}


string Collision::to_string(const Partition& p) const {
  stringstream ss;
  if (type == CollisionType::VOLMOL_VOLMOL) {
    ss << "coll_idx: " << colliding_molecule_id;
  }
  else {
    ss << "wall side: " << p.get_wall(colliding_wall_index).side;
  }

  ss << ", time: " << time << ", pos: " << pos;
  return ss.str();
}


void Collision::dump_array(Partition& p, const collision_vector_t& vec) {
  // printed in reverse - same as
  for (size_t i = 0; i < vec.size(); i++) {
    const char* str_type = (vec[i].type == CollisionType::VOLMOL_VOLMOL) ? "mol collision " : "wall collision ";
    cout << "  " << str_type << i << ": " << vec[i].to_string(p) << "\n";
  }
}

} /* namespace MCell */
