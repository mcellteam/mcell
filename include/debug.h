/******************************************************************************
 *
 * Copyright (C) 2020 by
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
