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

#include "molecule.h"
#include "world.h"

using namespace std;

namespace mcell {

// TODO: same as in dump_state.cpp, remove one of the copies
string get_molecule_flags_string(short flags) {
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
#undef DUMP_FLAG
	return res;
}

void base_molecule_t::dump_base(const std::string ind) const {
  cout << ind <<"flags: \t\t" << flags << "[uint16_t]\n";
  cout << ind <<"species_id: \t\t" << species_id << " [species_id_t]\n";
}

void volume_molecule_t::dump(const std::string ind) const {
  cout << ind <<"pos: \t\t" << pos << "[vec3_t]\n";
  cout << ind <<"subpartition_index: \t\t" << subpartition_index << " [uint32_t]\n";
  dump_base(ind);
}

void volume_molecule_t::dump(
		world_t* world,
		const std::string extra_comment,
		const std::string ind,
		const uint64_t iteration
) const {
	cout << ind << extra_comment << "it:" << iteration << ", idx:" << idx
			<< ", species " << world->species[species_id].name << ", pos:" << pos
			<< ", flags:" << get_molecule_flags_string(flags) << "\n";
}


} /* namespace mcell */
