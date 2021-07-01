/******************************************************************************
 *
 * Copyright (C) 2019-2021 by
 * The Salk Institute for Biological Studies
 *
 * Use of this source code is governed by an MIT-style
 * license that can be found in the LICENSE file or at
 * https://opensource.org/licenses/MIT.
 *
******************************************************************************/

#ifndef SRC4_SIMULATION_STATS_H_
#define SRC4_SIMULATION_STATS_H_

#include "defines.h"

namespace BNG {
class RxnContainer;
class RxnClass;
}

// for reporting of skipped reactions
struct RxnCountStats {
  RxnCountStats() :
    occurred(0),
    skipped(0) {
  }

  long long occurred;
  double skipped;
};

typedef std::map<std::string, RxnCountStats> NameStatsMap;
typedef std::map<std::string, NameStatsMap> BimolRxnCountStatsMap;

namespace MCell {

/*
 * Stats collected during simulation, contains also the number of the current iteration
 */
class SimulationStats {
public:
  SimulationStats() {
    reset(true);
  }
  void inc_ray_voxel_tests() {
    ray_voxel_tests++;
  }
  void inc_ray_polygon_tests() {
    ray_polygon_tests++;
  }
  void inc_ray_polygon_colls() {
    ray_polygon_colls++;
  }

  void inc_mol_wall_reflections() {
    mol_wall_reflections++;
  }

  void inc_vol_mol_vol_mol_collisions() {
    vol_mol_vol_mol_collisions++;
  }

  // new mcell4 stats
  void inc_mol_moves_between_walls() {
    mol_moves_between_walls++;
  }

  void inc_recomputations_of_counted_volume() {
    recomputations_of_counted_volume++;
  }

  void inc_num_waypoints_used() {
    num_waypoints_used++;
  }
  void inc_diffuse_3d_calls() {
    diffuse_3d_calls++;
  }

  void inc_diffusion_cummtime(const double steps) {
    diffusion_number++;
    diffusion_cummtime += steps; // this is a bit weird, steps are not time
  }

  void print_report();

  const uint64_t& get_current_iteration() const {
    return current_iteration;
  }

  uint64_t& get_current_iteration() {
    return current_iteration;
  }

  // to be used only in initialization
  void set_current_iteration(const uint64_t it) {
    current_iteration = it;
  }

  void inc_rxn_occurred(
      const BNG::RxnContainer& all_rxns,
      const BNG::RxnClass* rxn_class,
      const uint64_t occurred = 1);

  void inc_rxn_skipped(
      const BNG::RxnContainer& all_rxns,
      const BNG::RxnClass* rxn_class,
      const double skipped);

  void reset(const bool reset_also_current_iteration) {
    if (reset_also_current_iteration) {
      current_iteration = 0;
    }
    ray_voxel_tests = 0;
    ray_polygon_tests = 0;
    ray_polygon_colls = 0;
    mol_wall_reflections = 0;
    vol_mol_vol_mol_collisions = 0;
    mol_moves_between_walls = 0;
    num_waypoints_used = 0;
    recomputations_of_counted_volume = 0;
    diffuse_3d_calls = 0;
    diffusion_number = 0;
    diffusion_cummtime = 0.0;
  }

private:
  RxnCountStats& get_rxn_stats(
      const BNG::RxnContainer& all_rxns, const BNG::RxnClass* rxn_class);

  void print_missed_rxns_warnings();

  uint64_t current_iteration;

  // the stats below are not stored into checkpoint
  uint64_t ray_voxel_tests;
  uint64_t ray_polygon_tests;
  uint64_t ray_polygon_colls;

  uint64_t mol_wall_reflections;
  uint64_t vol_mol_vol_mol_collisions;

  uint64_t mol_moves_between_walls;

  uint64_t num_waypoints_used;
  uint64_t recomputations_of_counted_volume;

  uint64_t diffuse_3d_calls;

  uint64_t diffusion_number;
  double diffusion_cummtime;

  BimolRxnCountStatsMap bimol_rxn_stats;
};

} // namespace MCell

#endif // SRC4_SIMULATION_CONFIG_H_
