/******************************************************************************
 *
 * Copyright (C) 2019,2020 by
 * The Salk Institute for Biological Studies
 *
 * Use of this source code is governed by an MIT-style
 * license that can be found in the LICENSE file or at
 * https://opensource.org/licenses/MIT.
 *
******************************************************************************/

#ifndef SRC4_SIMULATION_CONFIG_H_
#define SRC4_SIMULATION_CONFIG_H_

#include "defines.h"

const char* const RUN_REPORT_PREFIX = "run_report_";

namespace MCell {

namespace API {
enum class WarningLevel;
};

/*
 * Constant data set in initialization useful for all classes, single object is owned by world
 */
class SimulationConfig: public BNG::BNGConfig {
public:
  SimulationConfig() :
    initial_time(TIME_INVALID),
    initial_iteration(UINT_INVALID),
    vacancy_search_dist2(FLT_INVALID),
    partition_edge_length(FLT_INVALID),
    num_subparts_per_partition_edge(UINT_INVALID),
    num_subparts_per_partition_edge_squared(UINT_INVALID),
    num_subparts(UINT_INVALID),
    subpart_edge_length(FLT_INVALID),
    subpart_edge_length_rcp(FLT_INVALID),
    num_radial_subdivisions(1024),
    use_expanded_list(true),
    randomize_smol_pos(false),
    check_overlapped_walls(true),
    rxn_class_cleanup_periodicity(0),
    species_cleanup_periodicity(0),
    molecules_order_random_shuffle_periodicity(DEFAULT_MOL_ORDER_SHUFFLE_PERIODICITY),
    sort_mols_by_subpart(false),
    memory_limit_gb(-1),
    iteration_report(true),
    wall_overlap_report(false),
    simulation_stats_every_n_iterations(0),
    continue_after_sigalrm(false),
    has_intersecting_counted_objects(false)
  {
    // enable debug assertions in libBNG
    debug_requires_diffusion_constants = true;
  }

  // configuration
  double initial_time; // simulation start time in us, non-zero if starting from a checkpoint
  uint64_t initial_iteration; // initial iteration, non-zero if starting from a checkpoint

  pos_t vacancy_search_dist2; /* Square of distance to search for free grid
                                  location to place surface product */

  Vec3 partition0_llf;

  pos_t partition_edge_length;
  uint num_subparts_per_partition_edge;
  uint num_subparts_per_partition_edge_squared;
  uint num_subparts; // == num_subpartitions_per_partition_edge^3
  pos_t subpart_edge_length; // == partition_edge_length / subpartitions_per_partition_dimension
  pos_t subpart_edge_length_rcp; // == 1/subpartition_edge_length

  uint num_radial_subdivisions; /* Size of 3D step length lookup tables, not configurable by user yet */
  std::vector<double> radial_2d_step; /* Lookup table of 2D diffusion step lengths (r_step_surface) */
  std::vector<double> radial_3d_step; /* Lookup table of 3D diffusion step lengths (r_step) */

  // other options
  bool use_expanded_list; /* If set, check neighboring subvolumes for mol-mol
                            interactions */
  bool randomize_smol_pos; /* If set, always place surface molecule at random
                             location instead of center of grid */
  bool check_overlapped_walls; /* Check geometry for overlapped walls? */

  API::WarningLevel molecule_placement_failure;

  uint rxn_class_cleanup_periodicity;
  uint species_cleanup_periodicity;
  uint molecules_order_random_shuffle_periodicity;

  bool sort_mols_by_subpart;

  int memory_limit_gb; // -1 means that limit is disabled

  // similar to MCell3's ITERATION_REPORT
  bool iteration_report;

  bool wall_overlap_report;

  int simulation_stats_every_n_iterations;

  bool continue_after_sigalrm;

  // initialized in World::init_counted_volumes
  // also tells whether waypoints in a partition were initialized
  bool has_intersecting_counted_objects;

  void init() {
    BNGConfig::init();
    init_subpartition_edge_length();
    init_radial_steps();
  }

  void dump();

  double get_simulation_start_time() const {
    assert(initial_time != TIME_INVALID);
    // simulation starts always in integer values of internal time
    double res = floor_to_multiple_f(initial_time / time_unit, 1);
    assert((int)res == res);
    return res;
  }

  // returns 'checkpoints/seed_<SEED>/it_' - without the iteration number
  std::string get_default_checkpoint_dir_prefix() const;

  std::string get_run_report_file_name() const;
  void initialize_run_report_file();

private:
  void init_subpartition_edge_length();
  void init_radial_steps();

};

std::string get_seed_dir_name(const int seed);

} // namespace MCell

#endif // SRC4_SIMULATION_CONFIG_H_
