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

class world_t;

#if 0
#define TYPE_SURF 0x001
#define TYPE_VOL 0x002
#define TYPE_MASK 0x003

#define ACT_DIFFUSE 0x008
#define ACT_REACT 0x020
#define ACT_NEWBIE 0x040
#define ACT_CHANGE 0x080
#define ACT_CLAMPED 0x1000

/* Flags telling us which linked lists the molecule appears in. */
#define IN_SCHEDULE 0x100
#define IN_SURFACE 0x200
#define IN_VOLUME 0x400
/* And a mask to pick off all three IN_ flags */
#define IN_MASK 0x700

/* Flags telling us what our counting status is */
#define COUNT_ME 0x800

/* Flag indicating that a molecule is old enough to take the maximum timestep */
#define MATURE_MOLECULE 0x2000
#endif

enum molecule_flags_e {
	MOLECULE_FLAG_DEFUNCT = 1 << 15,
};

class base_molecule_t {
public:
	base_molecule_t(const molecule_idx_t idx_, const species_id_t species_id_)
		: idx(idx_), flags(0), species_id(species_id_) {
	}

	molecule_idx_t idx; // necessary, index to partition's volume_molecules array
	uint16_t flags; // defunct by default, use bitfield instead?
	species_id_t species_id;

	bool is_defunct() const {
		return flags & MOLECULE_FLAG_DEFUNCT;
	}

	void set_is_defunct() {
		assert(!is_defunct() && "We really should not be defuncting one molecule multiple times");
		flags |= MOLECULE_FLAG_DEFUNCT;
	}

	// not using virtual methods, we do not want virtual methods table to be created for this object
	void dump_base(const std::string ind) const;
};


class volume_molecule_t : public base_molecule_t {
public:
	volume_molecule_t(const molecule_idx_t idx_, const species_id_t species_id_, const vec3_t& pos_)
		: base_molecule_t(idx_, species_id_),
			pos(pos_),
			subpartition_index(SUBPARTITION_INDEX_INVALID) {
	}

	vec3_t pos;
	uint32_t subpartition_index;

	void dump(const std::string ind) const;
	void dump(
			world_t* world,
			const std::string extra_comment,
			const std::string ind,
			const uint64_t iteration
	) const;

};

} // namespace mcell

#endif /* SRC4_MOLECULE_H_ */
