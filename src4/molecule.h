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

#ifndef SRC4_MOLECULE_H_
#define SRC4_MOLECULE_H_

#include "defines.h"

namespace mcell {


enum molecule_flags_e {
	MOLECULE_FLAG_DEFUNCT = 1 << 0,
};

class base_molecule_t {
public:
	base_molecule_t(const species_id_t species_id_)
		: flags(0), species_id(species_id_) {
	}

	uint16_t flags; // defunct by default, use bitfield instead?
	species_id_t species_id;

	bool is_defunct() {
		return flags & MOLECULE_FLAG_DEFUNCT;
	}

};


class volume_molecule_t : public base_molecule_t {
public:
	volume_molecule_t(const species_id_t species_id_, const vec3_t& pos_)
		: base_molecule_t(species_id_),
			pos(pos_),
			subpartition_index(SUBPARTITION_INDEX_INVALID) {
	}

	vec3_t pos;
	uint32_t subpartition_index;
};

} // namespace mcell

#endif /* SRC4_MOLECULE_H_ */
