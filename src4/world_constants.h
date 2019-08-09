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

#ifndef SRC4_WORLD_CONSTANTS_H_
#define SRC4_WORLD_CONSTANTS_H_

#include "defines.h"
#include "species.h"

namespace mcell {

/*
 * Constant data set in initialization useful for all classes, single object is owned by world
 */
struct world_constants_t {
  // configuration
  float_t time_unit;
  float_t length_unit;
  float_t rx_radius_3d;
  float_t partition_edge_length;
  uint subpartitions_per_partition_dimension;
  uint subpartitions_per_partition_dimension_squared;
  float_t subpartition_edge_length; // == partition_edge_length / subpartitions_per_partition_dimension
  float_t subpartition_edge_length_rcp; // == 1/subpartition_edge_length

  // other options
  bool use_expanded_list;


  const unimolecular_reactions_map_t* unimolecular_reactions_map; // owned by world
  const bimolecular_reactions_map_t* bimolecular_reactions_map; // owned by world

private:
  const std::vector<species_t>* species; // owned by world

private:
  void init_subpartition_edge_length() {
    if (partition_edge_length != 0) {
      subpartition_edge_length = partition_edge_length / (float_t)subpartitions_per_partition_dimension;
      subpartition_edge_length_rcp = 1.0/subpartition_edge_length;
    }
    subpartitions_per_partition_dimension_squared = powu(subpartitions_per_partition_dimension, 2);
  }

public:
  // called from world::init_simulation()
  void init(
      unimolecular_reactions_map_t* unimolecular_reactions_map_,
      bimolecular_reactions_map_t* bimolecular_reactions_map_,
      const std::vector<species_t>* species_
      ) {
    unimolecular_reactions_map = unimolecular_reactions_map_;
    bimolecular_reactions_map = bimolecular_reactions_map_;
    species = species_;
    init_subpartition_edge_length();
  }

  const species_t& get_species(species_id_t id) const {
    assert(id < species->size());
    return (*species)[id];
  }

  void dump();

  // TODO_LATER1: maybe add: bool fully_initialized;
};

} // namespace mcell

#endif // SRC4_WORLD_CONSTANTS_H_
