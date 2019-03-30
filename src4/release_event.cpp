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
extern "C" {
#include "rng.h" // MCell 3
#include "isaac64.h"
#include "mcell_structs.h"
}

#include <iostream>

#include "release_event.h"
#include "world.h"
#include "partition.h"

using namespace std;

namespace mcell {

void release_event_t::dump(const string ind) {
  cout << "Release event:\n";
  string ind2 = ind + "  ";
  base_event_t::dump(ind2);
  cout << ind2 << "location: \t\t" << location << " [vec3_t] \t\t\n";
  cout << ind2 << "species_id: \t\t" << species_id << " [species_id_t] \t\t\n";
  cout << ind2 << "release_number: \t\t" << release_number << " [uint32_t] \t\t\n";
  cout << ind2 << "name: \t\t" << name << " [string] \t\t\n";
}


void release_event_t::step() {
  // for now, let's simply release 'release_number' of molecules of 'species_id'
  // at 'location'
  partition_t& p = world->partitions[world->get_or_add_partition_index(location)];
  float_t time_step = world->species[species_id].time_step;
  uint32_t time_step_index = p.get_or_add_molecule_list_index_for_time_step(time_step);

  const int is_spheroidal = (release_shape == SHAPE_SPHERICAL ||
                             release_shape == SHAPE_ELLIPTIC ||
                             release_shape == SHAPE_SPHERICAL_SHELL);

  for (uint32_t i = 0; i < release_number; i++) {
    // TODO: this might use some refactoring
    vec3_t pos;
    do /* Pick values in unit square, toss if not in unit circle */
    {
      pos.x = (rng_dbl(&world->rng) - 0.5);
      pos.y = (rng_dbl(&world->rng) - 0.5);
      pos.z = (rng_dbl(&world->rng) - 0.5);
    } while (is_spheroidal &&
             pos.x * pos.x + pos.y * pos.y + pos.z * pos.z >= 0.25);

    if (release_shape == SHAPE_SPHERICAL_SHELL) {
      float_t r = sqrt(pos.x * pos.x + pos.y * pos.y + pos.z * pos.z) * 2.0;
      if (r == 0.0) {
        pos.x = 0.0;
        pos.y = 0.0;
        pos.z = 0.5;
      } else {
        pos /= r;
      }
    }

    float_t base_location[1][4];
    base_location[0][0] = pos.x * diameter.x + location.x;
    base_location[0][1] = pos.y * diameter.y + location.y;
    base_location[0][2] = pos.z * diameter.z + location.z;
    base_location[0][3] = 1;

    //TODO: t_matrix can be only identyty matrix for now, also use glm matrix mult.
    // mult_matrix(location, req->t_matrix, location, 1, 4, 4);

    vec3_t molecule_location;
    molecule_location.x = base_location[0][0];
    molecule_location.y = base_location[0][1];
    molecule_location.z = base_location[0][2];

    volume_molecule_t& new_vm = p.add_volume_molecule_with_time_step_index(
        volume_molecule_t(MOLECULE_ID_INVALID, species_id, molecule_location), time_step_index
    );
    new_vm.flags = TYPE_VOL | IN_VOLUME | ACT_DIFFUSE;
  }
}

} // namespace mcell


