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

#ifndef SRC4_SIMULATION_CONFIG_H_
#define SRC4_SIMULATION_CONFIG_H_

#include "defines.h"

namespace MCell {

/*
 * Constant data set in initialization useful for all classes, single object is owned by world
 */
// TODO: cleanup all unnecessary argument passing, e.g. in diffuse_react_event.cpp
class SimulationConfig: public BNG::BNGConfig {
public:
  SimulationConfig() :
    vacancy_search_dist2(FLT_INVALID),
    partition_edge_length(FLT_INVALID),
    num_subpartitions_per_partition(UINT_INVALID),
    num_subpartitions_per_partition_squared(UINT_INVALID),
    subpartition_edge_length(FLT_INVALID),
    subpartition_edge_length_rcp(FLT_INVALID),
    num_radial_subdivisions(1024),
    use_expanded_list(true),
    randomize_smol_pos(false),
    check_overlapped_walls(true),
    sort_mols_by_subpart(false),
    has_intersecting_counted_objects(false) {
  }

  // configuration

  float_t vacancy_search_dist2; /* Square of distance to search for free grid
                                  location to place surface product */

  Vec3 partition0_llf;

  float_t partition_edge_length; // TODO: rename to side
  uint num_subpartitions_per_partition; // TODO: rename to subpart...
  uint num_subpartitions_per_partition_squared;
  float_t subpartition_edge_length; // == partition_edge_length / subpartitions_per_partition_dimension
  float_t subpartition_edge_length_rcp; // == 1/subpartition_edge_length

  uint num_radial_subdivisions; /* Size of 3D step length lookup tables, not configurable by user yet */
  std::vector<float_t> radial_2d_step; /* Lookup table of 2D diffusion step lengths (r_step_surface) */
  std::vector<float_t> radial_3d_step; /* Lookup table of 3D diffusion step lengths (r_step) */

  // other options
  bool use_expanded_list; /* If set, check neighboring subvolumes for mol-mol
                            interactions */
  bool randomize_smol_pos; /* If set, always place surface molecule at random
                             location instead of center of grid */
  bool check_overlapped_walls; /* Check geometry for overlapped walls? */

  bool sort_mols_by_subpart;

  // initialized in World::init_counted_volumes
  // also tells whether waypoints in a partition were initialized
  bool has_intersecting_counted_objects;


  void init() {
    BNGConfig::init();
    init_subpartition_edge_length();
    init_radial_steps();
  }

  void dump();

private:
  void init_subpartition_edge_length();
  void init_radial_steps();

};

} // namespace MCell

#endif // SRC4_SIMULATION_CONFIG_H_
