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

#include "mcell_structs.h"

#include "molecule.h"
#include "geometry.h"
#include "world.h"

using namespace std;


namespace MCell {

// TODO_LATER: same as in dump_state.cpp, remove one of the copies
// for now, dump_state might be useful also without mcell4
string get_molecule_flags_string(uint flags, bool full_dump = true) {
  string res;
#define DUMP_FLAG(f, mask) if (((f) & (mask)) != 0) res += string(#mask) + ", ";
  DUMP_FLAG(flags, TYPE_SURF)
  DUMP_FLAG(flags, TYPE_VOL)
  DUMP_FLAG(flags, ACT_DIFFUSE)
  if (full_dump) {
    DUMP_FLAG(flags, ACT_REACT)

    if ((flags & MOLECULE_FLAG_RESCHEDULE_UNIMOL_RX) != 0)
      res += string("ACT_NEWBIE") + ", ";

    if ((flags & MOLECULE_FLAG_RESCHEDULE_UNIMOL_RX) != 0)
      res += string("ACT_CHANGE") + ", ";
  }
  DUMP_FLAG(flags, ACT_CLAMPED)
  if (full_dump) {
    DUMP_FLAG(flags, IN_SCHEDULE)
    DUMP_FLAG(flags, IN_SURFACE)
  }
  DUMP_FLAG(flags, IN_VOLUME)
  DUMP_FLAG(flags, MOLECULE_FLAG_DEFUNCT)
#undef DUMP_FLAG
  return res;
}


void Molecule::dump(const string ind) const {
  if (is_vol()) {
    cout << ind << "pos: \t\t" << v.pos << " [vec3_t]\n";
    cout << ind << "subpartition_index: \t\t" << v.subpart_index << " [uint]\n";
  }
  else if (is_surf()) {
    cout << ind << "pos: \t\t" << s.pos << " [vec2_t]\n";
  }
  cout << ind << "flags: \t\t" << flags << " [uint]\n";
  cout << ind << "species_id: \t\t" << species_id << " [species_id_t]\n";
}


void Molecule::dump(
    const Partition& p,
    const string extra_comment,
    const string ind,
    const uint64_t iteration,
    const float_t time,
    const bool print_position
) const {
  cout
    << ind << extra_comment << "it:" << iteration << ", idx:" << id
    << ", species: " << p.get_world_constants().get_species(species_id).name;

  if (print_position) {
    cout << ", pos:";

    if (is_vol()) {
      cout << v.pos;
    }
    else if (is_surf()) {
      cout << s.pos;
      const Wall& w = p.get_wall(s.wall_index);
      cout << ", wall side: " << w.side;
      cout << ", grid index: " << s.grid_tile_index;
    }
  }
  cout
    << ", flags:" << get_molecule_flags_string(flags, false);
#ifdef DEBUG_SUBPARTITIONS
  if (is_vol()) {
    cout << ", subpartition:" << subpart_index;
  }
#endif

  cout
      << ", time: " << time
      << "\n";
}


string Molecule::to_string() const {
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
  ss << ", flags:" << get_molecule_flags_string(flags, false);
  return ss.str();
}


void Molecule::dump_array(const std::vector<Molecule>& vec) {
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
