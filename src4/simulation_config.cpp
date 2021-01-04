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
#include <cmath>

using namespace std;

namespace MCell {

void SimulationConfig::init_subpartition_edge_length() {
  release_assert(partition_edge_length > 0);
  subpartition_edge_length = partition_edge_length / (float_t)num_subpartitions_per_partition_edge;
  subpartition_edge_length_rcp = 1.0/subpartition_edge_length;
  num_subpartitions_per_partition_edge_squared = powu(num_subpartitions_per_partition_edge, 2);
  num_subpartitions = powu(num_subpartitions_per_partition_edge, 3);
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
/***************************************************************************
init_r_step_surface:
  In: number of desired radial subdivisions
  Out: pointer to array of doubles containing those subdivisions
       returns NULL on malloc failure
  Note: This is for 3D molecules emitted from a plane
***************************************************************************/
void SimulationConfig::init_radial_steps() {
  float_t inc, target, accum, r, r_max, delta_r, delta_r2;

  // 3D
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

  // 2D
  radial_2d_step.resize(num_radial_subdivisions);
  static const float_t sqrt_pi = 1.7724538509055160273;

  float_t step = 1.0 / num_radial_subdivisions;
  int i = 0;
  float_t p = (1.0 - 1e-6) * step;
  r = 0;
  for (; p < 1.0; p += step, i++) {
    float_t r_min = 0;
    float_t r_max = 3.0;          /* 17 bit high-end CDF cutoff */
    for (int j = 0; j < 20; j++) /* 20 bits of accuracy */
    {
      r = 0.5 * (r_min + r_max);
      float_t cdf = 1.0 - exp(-r * r) + sqrt_pi * r * erfc(r);
      if (cdf > p)
        r_max = r;
      else
        r_min = r;
    }
    radial_2d_step[i] = r;
  }
}



void SimulationConfig::dump() {
  BNGConfig::dump();
  cout << "SimulationConfig:\n";
#define DUMP_ATTR(A) cout << "  " #A ": \t\t" << A << "\n"
  DUMP_ATTR(vacancy_search_dist2);
  DUMP_ATTR(partition0_llf);
  DUMP_ATTR(partition_edge_length);
  DUMP_ATTR(num_subpartitions_per_partition_edge);
  DUMP_ATTR(num_subpartitions_per_partition_edge_squared);
  DUMP_ATTR(num_subpartitions);
  DUMP_ATTR(subpartition_edge_length);
  DUMP_ATTR(subpartition_edge_length_rcp);
  DUMP_ATTR(num_radial_subdivisions);
  DUMP_ATTR(use_expanded_list);
  DUMP_ATTR(randomize_smol_pos);
  DUMP_ATTR(check_overlapped_walls);
  DUMP_ATTR(rxn_class_cleanup_periodicity);
  DUMP_ATTR(species_cleanup_periodicity);
  DUMP_ATTR(sort_mols_by_subpart);
  DUMP_ATTR(memory_limit_gb);
  DUMP_ATTR(simulation_stats_every_n_iterations);
  DUMP_ATTR(has_intersecting_counted_objects);
#undef DUMP_ATTR
}

} // namespace MCell
