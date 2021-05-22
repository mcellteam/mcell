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

#include "bng/species.h"

#include "mcell_structs_shared.h"

#include "molecule.h"
#include "geometry.h"
#include "world.h"

#include "debug_config.h"

using namespace std;


namespace MCell {

void Molecule::set_no_need_to_schedule_flag(
    const BNG::SpeciesContainer& all_species) {
  clear_flag(MOLECULE_FLAG_NO_NEED_TO_SCHEDULE);

  const BNG::Species& s = all_species.get(species_id);
  if (s.D != 0 ||
      s.has_flag(BNG::SPECIES_FLAG_HAS_UNIMOL_RXN) ||
      s.has_flag(BNG::SPECIES_FLAG_CAN_SURFSURF) ||
      s.has_flag(BNG::SPECIES_FLAG_CAN_SURFWALL) ||
      s.has_flag(BNG::SPECIES_FLAG_CAN_INTERMEMBRANE_SURFSURF)
  ) {
    return;
  }

  set_flag(MOLECULE_FLAG_NO_NEED_TO_SCHEDULE);
}


string get_molecule_flags_string(uint flags, bool full_dump = true) {
  string res;
#define DUMP_FLAG(f, mask) if (((f) & (mask)) != 0) res += string(#mask) + ", ";
  DUMP_FLAG(flags, MOLECULE_FLAG_SURF)
  DUMP_FLAG(flags, MOLECULE_FLAG_VOL)
  if (full_dump) {
    if ((flags & MOLECULE_FLAG_SCHEDULE_UNIMOL_RXN) != 0)
      res += string("ACT_NEWBIE") + ", ";

    if ((flags & MOLECULE_FLAG_SCHEDULE_UNIMOL_RXN) != 0)
      res += string("ACT_CHANGE") + ", ";
  }
  DUMP_FLAG(flags, MOLECULE_FLAG_ACT_CLAMPED)
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
    const double time,
    const bool print_position,
    const bool print_flags
) const {
  cout
    << ind << extra_comment << "it: " << iteration << ", id: " << id
    << ", species: " << p.get_species(species_id).name;

  if (print_position) {
    cout << ", pos: ";

    if (is_vol()) {
      cout << v.pos;
    }
    else if (is_surf()) {
      cout << s.pos;
      cout << ", orient: " << s.orientation;
      const Wall& w = p.get_wall(s.wall_index);
      cout << ", wall side: " << w.side;
      cout << ", grid index: " << s.grid_tile_index;
    }
  }
  if (print_flags) {
    cout  << ", flags: " << get_molecule_flags_string(flags, false);
  }

#ifdef DEBUG_SUBPARTITIONS
  IVec3 indices;
  p.get_subpart_3d_indices_from_index(v.subpart_index, indices);
  if (is_vol()) {
    cout << ", subpartition: " << v.subpart_index << " " << indices;
  }
#endif

#ifdef DEBUG_COUNTED_VOLUMES
  if (is_vol()) {
    cout << ", counted vols id: " << v.counted_volume_index;
    p.get_counted_volume(v.counted_volume_index).dump();
  }
#endif

#ifdef DEBUG_COMPARTMENTS
  cout << ", compartment id: " << reactant_compartment_id;
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
