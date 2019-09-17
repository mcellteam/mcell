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
#include "logging.h"
}

#include <iostream>

#include "release_event.h"
#include "world.h"
#include "partition.h"

#include "geometry_utils.inc"
#include "grid_utils.inc"

using namespace std;

namespace MCell {

void ReleaseEvent::dump(const string ind) {
  cout << "Release event:\n";
  string ind2 = ind + "  ";
  BaseEvent::dump(ind2);
  cout << ind2 << "location: \t\t" << location << " [vec3_t] \t\t\n";
  cout << ind2 << "species_id: \t\t" << species_id << " [species_id_t] \t\t\n";
  cout << ind2 << "release_number: \t\t" << release_number << " [uint] \t\t\n";
  cout << ind2 << "name: \t\t" << release_site_name << " [string] \t\t\n";
}


uint ReleaseEvent::calculate_number_to_release() {
  // release_number_method - only 0 allowed now
  // case CONSTNUM:
  return release_number;
}


// NOTE: maybe a template will be needed for this function, used a lot in mcell3
static size_t cum_area_bisect_high(const vector<CummAreaPWallIndexPair>& array, float_t val) {
  size_t low = 0;
  size_t hi = array.size() - 1;
  size_t mid = 0;

  while (hi - low > 1) {
    mid = (hi + low) / 2;
    if (array[mid].first > val) {
      hi = mid;
    } else {
      low = mid;
    }
  }

  if (array[low].first > val)
  {
    return low;
  }
  else {
    return hi;
  }
}


void ReleaseEvent::place_single_molecule_onto_grid(Partition& p, Wall& wall, tile_index_t tile_index) {

  vec2_t pos_on_wall = GridUtil::grid2uv_random(wall, tile_index, world->rng);

  Molecule& new_sm = p.add_surface_molecule(
      Molecule(MOLECULE_ID_INVALID, species_id, pos_on_wall)
  );

  new_sm.s.wall_index = wall.index;
  new_sm.s.orientation = orientation;

  new_sm.s.grid_tile_index = tile_index;
  wall.grid.set_molecule_tile(tile_index, new_sm.id);

  new_sm.flags = ACT_NEWBIE | TYPE_SURF | ACT_DIFFUSE | IN_SURFACE;
}


void ReleaseEvent::release_onto_regions(uint computed_release_number) {
  int success = 0, failure = 0;
  float_t seek_cost = 0;

  // TODO_LATER: for now we are assuming that we have just a single partition
  // and releases do not cross partition boundary

  assert(!cum_area_and_pwall_index_pairs.empty());
  float_t total_area = cum_area_and_pwall_index_pairs.back().first;
  float_t est_sites_avail = (int)total_area;
  const float_t rel_list_gen_cost = 10.0; /* Just a guess */
  float_t pick_cost = rel_list_gen_cost * est_sites_avail;

  uint n = computed_release_number;

  const int too_many_failures = 10; /* Just a guess */
  while (n > 0) {
    if (failure >= success + too_many_failures) {
      seek_cost =
          n * (((double)(success + failure + 2)) / ((double)(success + 1)));
    }

    if (seek_cost < pick_cost) {
      float_t A = rng_dbl(&world->rng) * total_area;
      size_t cum_area_index = cum_area_bisect_high(cum_area_and_pwall_index_pairs, A);
      PartitionWallIndexPair pw = cum_area_and_pwall_index_pairs[cum_area_index].second;
      Partition& p = world->get_partition(pw.first);
      Wall& wall = p.get_wall(pw.second);

      if (!wall.has_initialized_grid()) {
        wall.initialize_grid(p); // sets wall's grid_index
      }

      Grid& grid = wall.grid;

      // get the random number for the current wall
      if (cum_area_index != 0) {
        A -= cum_area_and_pwall_index_pairs[cum_area_index - 1].first;
      }

      tile_index_t tile_index = (grid.num_tiles_along_axis * grid.num_tiles_along_axis) * (A / wall.area);
      if (tile_index >= grid.num_tiles) {
        tile_index = grid.num_tiles - 1;
      }

      if (grid.get_molecule_on_tile(tile_index) != MOLECULE_ID_INVALID) {
        failure++;
        continue;
      }

      place_single_molecule_onto_grid(p, wall, tile_index);

      success++;
      n--;
    }
    else {
      assert(false && "Recovery from too many failures during surf mol release is not implemented yet");
    }
  }
}


void ReleaseEvent::release_ellipsoid_or_rectcuboid(uint computed_release_number) {

  Partition& p = world->get_partition(world->get_or_add_partition_index(location));
  float_t time_step = world->get_species(species_id).time_step;

  const int is_spheroidal = (release_shape == SHAPE_SPHERICAL ||
                             release_shape == SHAPE_ELLIPTIC ||
                             release_shape == SHAPE_SPHERICAL_SHELL);

  for (uint i = 0; i < computed_release_number; i++) {
    vec3_t pos;
    do /* Pick values in unit square, toss if not in unit circle */
    {
      pos.x = (rng_dbl(&world->rng) - 0.5);
      pos.y = (rng_dbl(&world->rng) - 0.5);
      pos.z = (rng_dbl(&world->rng) - 0.5);
    } while (is_spheroidal && len3_squared(pos) >= 0.25);

    if (release_shape == SHAPE_SPHERICAL_SHELL) {
      float_t r = sqrt(len3_squared(pos)) * 2.0;
      if (r == 0.0) {
        pos = vec3_t(0.0, 0.0, 0.5);
      } else {
        pos /= r;
      }
    }

    float_t base_location[1][4];
    base_location[0][0] = pos.x * diameter.x + location.x;
    base_location[0][1] = pos.y * diameter.y + location.y;
    base_location[0][2] = pos.z * diameter.z + location.z;
    base_location[0][3] = 1;

    // TODO_LATER: t_matrix can be only identity matrix for now, also use glm matrix mult.
    // mult_matrix(location, req->t_matrix, location, 1, 4, 4);

    vec3_t molecule_location;
    molecule_location.x = base_location[0][0];
    molecule_location.y = base_location[0][1];
    molecule_location.z = base_location[0][2];

    // TODO_LATER: location can be close to a partition boundary, we might need to release to a different partition
    Molecule& new_vm = p.add_volume_molecule(
        Molecule(MOLECULE_ID_INVALID, species_id, molecule_location)
    );
    new_vm.flags = ACT_NEWBIE | TYPE_VOL | IN_VOLUME | ACT_DIFFUSE;
  }
}


void ReleaseEvent::step() {
  // for now, let's simply release 'release_number' of molecules of 'species_id'
  // at 'location'

  uint number = calculate_number_to_release();

  const Species& species = world->get_species(species_id);

  if (release_shape == SHAPE_REGION) {
    if (species.is_surf()) {
      release_onto_regions(number);
    }
    else {
      assert(false && "Region volume mol release is not supported yet.");
    }
  }
  else {
    assert(diameter.is_valid());
    release_ellipsoid_or_rectcuboid(number);
  }

  cout
    << "Released " << number << " " << species.name << " from \"" << release_site_name << "\""
    << " at iteration " << world->current_iteration << ".\n";
}


} // namespace mcell


