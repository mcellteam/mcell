/******************************************************************************
 *
 * Copyright (C) 2019 by
 * The Salk Institute for Biological Studies and
 * Pittsburgh Supercomputing Center, Carnegie Mellon University
 *
 * Use of this source code is governed by an MIT-style
 * license that can be found in the LICENSE file or at
 * https://opensource.org/licenses/MIT.
 *
******************************************************************************/

#include <iostream>
#include <string>
#include <sstream>

#include "collision_structs.h"

#include "partition.h"
#include "molecule.h"

#include "debug_config.h"

using namespace std;

namespace MCell {


void Collision::dump(Partition& p, const std::string ind) const {
  cout << ind << "diffused_molecule:\n";
  p.get_m(diffused_molecule_id).dump(ind + "  ");
  if (type == CollisionType::VOLMOL_VOLMOL) { // dump is rather limited for now, goes not deal with all types
    cout << ind << "colliding_molecule:\n";
    p.get_m(colliding_molecule_id).dump(ind + "  ");
    cout << ind << "reaction:";

    if (rxn_class != nullptr) {
      rxn_class->dump(ind + "  ");
    }
  }
  else {
    cout << ind << "colliding_wall_index: " << colliding_wall_index << "\n";
  }

  cout << "time: \t\t" << time << " [double] \t\t\n";
  cout << "position: \t\t" << pos << " [vec3_t] \t\t\n";
}

void Collision::dump(
    const Partition& p,
    const std::string extra_comment,
    const uint64_t iteration,
    double time_override
) const {

  cout << extra_comment << "it:" << iteration << ", ";

  if (is_mol_mol_reaction()) {
    cout <<
      "bimol rxn" <<
      ", idA:"  << diffused_molecule_id <<
      ", idB:"  << colliding_molecule_id <<
      ", time: " << ((time_override == TIME_INVALID) ? time : time_override);

      if (type != CollisionType::SURFMOL_SURFMOL) {
        cout << ", pos " << pos;
      }
  }
  else if (is_unimol_reaction()) {
    cout <<
      "unimol rxn" <<
      ", idA:"  << diffused_molecule_id <<
      ", time: " << ((time_override == TIME_INVALID) ? time : time_override);
  }
  else if (is_wall_collision()) {
    #ifndef NODEBUG_WALL_COLLISIONS
      cout << "wall collision" <<
          ", idA:"  << diffused_molecule_id <<
          //", wall index:"  << colliding_wall_index <<
          ", time: " << ((time_override == TIME_INVALID) ? time : time_override) <<
          ", pos " << pos;
    #endif
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
    ss << "obj name: " << p.get_geometry_object(p.get_wall(colliding_wall_index).object_index).name;
    ss << ", wall index: " << colliding_wall_index;
    ss << ", wall side: " << p.get_wall(colliding_wall_index).side;
#ifdef DEBUG_COUNTED_VOLUMES
    ss << ", geom obj id: " << p.get_wall(colliding_wall_index).object_id;
#endif
  }

  ss << ", time: " << time << ", pos: " << pos;
  return ss.str();
}


void Collision::dump_array(Partition& p, const CollisionsVector& vec) {
  // printed in reverse - same as
  for (size_t i = 0; i < vec.size(); i++) {
    #ifdef NODEBUG_WALL_COLLISIONS
      if (vec[i].type == CollisionType::WALL_FRONT || vec[i].type == CollisionType::WALL_BACK) {
        continue;
      }
    #endif

    const char* str_type = (vec[i].type == CollisionType::VOLMOL_VOLMOL) ? "mol collision " : "wall collision ";
    cout << "  " << str_type << i << ": " << vec[i].to_string(p) << "\n";
  }
}

} /* namespace MCell */
