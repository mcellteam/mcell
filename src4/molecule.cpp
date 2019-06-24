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
#include <string>
#include <sstream>

#include "molecule.h"
#include "world.h"

using namespace std;


namespace mcell {

// TODO: same as in dump_state.cpp, remove one of the copies
static string get_molecule_flags_string(uint32_t flags) {
  string res;
#define DUMP_FLAG(f, mask) if (((f) & (mask)) != 0) res += string(#mask) + ", ";
  DUMP_FLAG(flags, TYPE_SURF)
  DUMP_FLAG(flags, TYPE_VOL)
  DUMP_FLAG(flags, ACT_DIFFUSE)
  DUMP_FLAG(flags, ACT_REACT)
  DUMP_FLAG(flags, ACT_NEWBIE)
  DUMP_FLAG(flags, ACT_CHANGE)
  DUMP_FLAG(flags, ACT_CLAMPED)
  DUMP_FLAG(flags, IN_SCHEDULE)
  DUMP_FLAG(flags, IN_SURFACE)
  DUMP_FLAG(flags, IN_VOLUME)
  DUMP_FLAG(flags, MOLECULE_FLAG_DEFUNCT)
#undef DUMP_FLAG
  return res;
}


void molecule_t::dump(const string ind) const {
  if (is_vol()) {
    cout << ind << "pos: \t\t" << v.pos << " [vec3_t]\n";
    cout << ind << "subpartition_index: \t\t" << v.subpart_index << " [uint32_t]\n";
  }
  else if (is_surf()) {
    cout << ind << "pos: \t\t" << s.pos << " [vec2_t]\n";
  }
  cout << ind << "flags: \t\t" << flags << " [uint32_t]\n";
  cout << ind << "species_id: \t\t" << species_id << " [species_id_t]\n";
}


void molecule_t::dump(
    const world_t* world,
    const string extra_comment,
    const string ind,
    const uint64_t iteration,
    const float_t time
) const {
  cout
    << ind << extra_comment << "it:" << iteration << ", idx:" << id
    << ", species " << world->species[species_id].name << ", pos:";

  if (is_vol()) {
    cout << v.pos;
  }
  else if (is_surf()) {
    cout << s.pos;
  }
  cout
    << ", flags:" << get_molecule_flags_string(flags);
#ifdef DEBUG_SUBPARTITIONS
  if (is_vol()) {
    cout << ", subpartition:" << subpart_index;
  }
#endif

  cout
      << ", time: " << time
      << "\n";
}


string molecule_t::to_string() const {
  stringstream ss;
  ss
    << "id: " << id
    << ", species: " << species_id
    << ", pos: ";
  if (is_vol()) {
    cout << v.pos;
  }
  else if (is_surf()) {
    cout << s.pos;
  }
  ss << ", flags:" << get_molecule_flags_string(flags);
  return ss.str();
}


void molecule_t::dump_array(const std::vector<molecule_t>& vec) {
  for (size_t i = 0; i < vec.size(); i++) {
    if (vec[i].is_vol()) {
      cout << "  vm ";
    }
    else if (vec[i].is_surf()) {
      cout << "  sm ";
    }
    cout << i << ": " << vec[i].to_string() << "\n";
  }
}

} // namespace mcell
