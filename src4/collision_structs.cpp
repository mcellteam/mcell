/*
 * collision_structs.cpp
 *
 *  Created on: Dec 4, 2019
 *      Author: ahusar
 */


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
  if (type == CollisionType::VOLMOL_VOLMOL) {
    cout << ind << "colliding_molecule:\n";
    p.get_m(colliding_molecule_id).dump(ind + "  ");
    cout << ind << "reaction:";
    rx->dump(ind + "  ");
  }
  else {
    cout << ind << "colliding_wall_index: " << colliding_wall_index << "\n";
  }

  cout << "time: \t\t" << time << " [float_t] \t\t\n";
  cout << "position: \t\t" << pos << " [vec3_t] \t\t\n";
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
