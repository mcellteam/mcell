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

#include <iostream>

#include "bng/species.h"
#include "bng/bng_defines.h"

#include <numeric>      // std::iota
#include <algorithm>    // std::sort, std::stable_sort

using namespace std;

namespace BNG {

// based on assemble_mol_species
void Species::update_space_and_time_step(const float_t time_unit, const float_t length_unit) {
  // Immobile (boring)
  if (!distinguishable_f(D, 0, EPS)) {
    space_step = 0.0;
    time_step = 1.0;
  }
  // Custom timestep or spacestep
  // not supported yet
  /*else if (new_spec->time_step != 0.0) {
    // Hack--negative value means custom space step
    if (new_spec->time_step < 0) {
      double lr_bar = -new_spec->time_step;
      if (species->is_2d) {
        new_spec->time_step =
            lr_bar * lr_bar / (MY_PI * 1.0e8 * new_spec->D * global_time_unit);
      } else {
        new_spec->time_step =
            lr_bar * lr_bar * MY_PI /
            (16.0 * 1.0e8 * new_spec->D * global_time_unit);
      }
      new_spec->space_step =
          sqrt(4.0 * 1.0e8 * new_spec->D * new_spec->time_step *
               global_time_unit) *
          state->r_length_unit;
    }
    else {
      new_spec->space_step =
          sqrt(4.0 * 1.0e8 * new_spec->D * new_spec->time_step) *
          state->r_length_unit;
      new_spec->time_step /= global_time_unit;
    }
  }*/
  // Global timestep (this is the typical case)
  else /*if (!distinguishable(state->space_step, 0, EPS_C))*/ {
    space_step = sqrt_f(4.0 * 1.0e8 * D * time_unit) / length_unit;
    time_step = 1.0;
  }
  /*// Global spacestep
  else {
    double space_step = state->space_step * state->length_unit;
    if (species->is_2d) {
      new_spec->time_step =
          space_step * space_step /
          (MY_PI * 1.0e8 * new_spec->D * global_time_unit);
    }
    else {
      new_spec->time_step =
          space_step * space_step * MY_PI /
          (16.0 * 1.0e8 * new_spec->D * global_time_unit);
    }
    new_spec->space_step = sqrt(4.0 * 1.0e8 * new_spec->D *
                                new_spec->time_step * global_time_unit) *
                                state->r_length_unit;
  }*/
}


void Species::dump(const BNGData& bng_data, const string ind) const {
  cout << ind << "species_id: \t\t" << id << " [uint16_t] \t\t/* Unique ID for this species */\n";
  cout << ind << "name: *\t\t" << name << " [string] \t\t/* Symbol table entry (name) */\n";
  cout << ind << "D: \t\t" << D << " [float_t] \t\t/* Diffusion constant */\n";
  cout << ind << "space_step: \t\t" << space_step << " [float_t] \t\t/* Characteristic step length */\n";
  cout << ind << "time_step: \t\t" << time_step << " [float_t] \t\t/* Minimum (maximum?) sensible timestep */\n";
  cout << ind << "flags: \t\t0x" << hex << get_flags() << dec << " [uint] \t\t/* Flags */\n";
  cout << ind << "CplxInstance: ";
  CplxInstance::dump(bng_data, true, ind + "  ");
  cout << "\n";
}


template <typename T>
vector<size_t> sort_indexes(const T &v) {

  // initialize original index locations
  vector<size_t> idx(v.size());
  iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  // using std::stable_sort instead of std::sort
  // to avoid unnecessary index re-orderings
  // when v contains elements of equal values
  stable_sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1].name < v[i2].name;});

  return idx;
}


void Species::dump_array(const BNGData& bng_data, const SpeciesVector& vec, const bool sorted) {
  cout << "Species array: " << (vec.empty() ? "EMPTY" : "") << "\n";

  if (sorted) {
    // dump sorted by name
    vector<size_t> sorted_indices = sort_indexes(vec);
    for (auto i: sorted_indices) {
      cout << i << ":\n";
      vec[i].dump(bng_data, "  ");
    }
  }
  else {
    for (size_t i = 0; i < vec.size(); i++) {
      cout << i << ":\n";
      vec[i].dump(bng_data, "  ");
    }
  }
}

} // namespace mcell
