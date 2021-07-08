/******************************************************************************
 *
 * Copyright (C) 2020 by
 * The Salk Institute for Biological Studies and
 * Pittsburgh Supercomputing Center, Carnegie Mellon University
 *
 * Use of this source code is governed by an MIT-style
 * license that can be found in the LICENSE file or at
 * https://opensource.org/licenses/MIT.
 *
******************************************************************************/

#include "defines.h"

namespace MCell {

// TODO: move to cpp
static void dump_vol_mol_timing(
    std::string extra_comment,
    uint64_t iteration,
    molecule_id_t id,
    double scheduled_time, double max_time, double unimol_time,
    double rate_factor, double r_rate_factor, double steps, double t_steps
) {

  std::cout
     << extra_comment << ": it:" << iteration << ", id:" << id
     << ", scheduled_time: " << scheduled_time
     << ", max_time: " << max_time
     << ", unimol_time: " << ((unimol_time != TIME_INVALID && unimol_time != scheduled_time) ? unimol_time : 0)
     << ", rate_factor: " << rate_factor
     << ", r_rate_factor: " << r_rate_factor
     << ", steps: " << steps
     << ", t_steps: " << t_steps
     << "\n";
}


static void dump_surf_mol_timing(
    std::string extra_comment,
    uint64_t iteration,
    molecule_id_t id,
    double scheduled_time, double max_time, double unimol_time,
    double space_factor, double steps, double t_steps
) {

  std::cout
     << extra_comment << ": it:" << iteration << ", id:" << id
     << ", scheduled_time: " << scheduled_time
     << ", max_time: " << max_time
     << ", unimol_time: " << ((unimol_time != TIME_INVALID && unimol_time != scheduled_time) ? unimol_time : 0)
     << ", space_factor: " << space_factor
     << ", steps: " << steps
     << ", t_steps: " << t_steps
     << "\n";
}


static void dump_react_2D_all_neighbors_timing(
    double time,
    double mol_time
) {
  std::cout
    << "react_2D_all_neighbors: "
    << "time: " << time
    << ", mol_time (sm->t): " << mol_time
    << "\n";
}

static void dump_outcome_bimolecular_timing(
    double time
) {
  std::cout
    << "outcome_bimolecular: time: " << time
    << "\n";
}

static void dump_uint_vector(const std::vector<uint> v) {
  for (uint i = 0; i < v.size(); i++) {
    std::cout << v[i] << ", ";
  }
  std::cout << "\n";
}

static void dump_uint_set(const std::set<uint> s) {
  for (uint val: s) {
    std::cout << val << ", ";
  }
  std::cout << "\n";
}

} // namespace?
