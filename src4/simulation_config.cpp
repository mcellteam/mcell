/******************************************************************************
 *
 * Copyright (C) 2019,2020 by
 * The Salk Institute for Biological Studies
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

#include "simulation_config.h"

#include <iostream>
using namespace std;

namespace MCell {

void SimulationConfig::init_subpartition_edge_length() {
  release_assert(partition_edge_length > 0);
  subpartition_edge_length = partition_edge_length / (float_t)num_subpartitions_per_partition;
  subpartition_edge_length_rcp = 1.0/subpartition_edge_length;
  num_subpartitions_per_partition_squared = powu(num_subpartitions_per_partition, 2);
}


/***************************************************************************
r_func:
  In: double containing distance (arbitrary units, mean=1.0)
  Out: double containing probability of diffusing that distance
***************************************************************************/
float_t r_func(float_t s) {
  float_t f, s_sqr, val;

  f = 2.2567583341910251478; /* 4.0/sqrt(pi) */
  s_sqr = s * s;
  val = f * s_sqr * exp(-s_sqr);

  return val;
}


/***************************************************************************
init_r_step:
  In: number of desired radial subdivisions
  Out: pointer to array of doubles containing those subdivisions
       returns NULL on malloc failure
  Note: This is for 3D diffusion from a point source (molecule movement)
***************************************************************************/
void SimulationConfig::init_radial_3d_step() {
  float_t inc, target, accum, r, r_max, delta_r, delta_r2;

  radial_3d_step.resize(num_radial_subdivisions);

  inc = 1.0 / num_radial_subdivisions;
  accum = 0;
  r_max = 3.5;
  delta_r = r_max / (1000 * num_radial_subdivisions);
  delta_r2 = 0.5 * delta_r;
  r = 0;
  target = 0.5 * inc;
  uint j = 0;
  while (j < num_radial_subdivisions) {
    accum = accum + (delta_r * r_func(r + delta_r2));
    r = r + delta_r;
    if (accum >= target) {
      radial_3d_step[j] = r;
      target = target + inc;
      j++;
    }
  }
}



void SimulationConfig::dump() {
  BNGConfig::dump();
  cout << "SimulationConfig:\n";
  cout << "  vacancy_search_dist2: \t\t" << vacancy_search_dist2 << " [float_t] \t\t\n";
  cout << "  partition_edge_length: \t\t" << partition_edge_length << " [float_t] \t\t\n";
  cout << "  num_subpartitions_per_partition: \t\t" << num_subpartitions_per_partition << " [uint] \t\t\n";
  cout << "  num_subpartitions_per_partition_squared: \t\t" << num_subpartitions_per_partition_squared << " [uint] \t\t\n";
  cout << "  subpartition_edge_length: \t\t" << subpartition_edge_length << " [float_t] \t\t\n";
  cout << "  subpartition_edge_length_rcp: \t\t" << subpartition_edge_length_rcp << " [float_t] \t\t\n";
  cout << "  use_expanded_list: \t\t" << use_expanded_list << " [bool] \t\t\n";
  cout << "  randomize_smol_pos: \t\t" << randomize_smol_pos << " [bool] \t\t\n";
}


} // namespace MCell
