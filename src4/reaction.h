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


#ifndef SRC4_REACTION_H_
#define SRC4_REACTION_H_

#include "defines.h"

namespace mcell {

typedef int32_t orientation_t;
const orientation_t ORIENTATION_DOWN = -1;
const orientation_t ORIENTATION_NONE = 0;
const orientation_t ORIENTATION_UP = 1;

struct species_with_orientation_t {
	species_with_orientation_t()
		: species_id(SPECIES_ID_INVALID), orientation(ORIENTATION_NONE) {
	}
	species_with_orientation_t(
			const species_id_t species_id_,
			const orientation_t orientation_)
		: species_id(species_id_), orientation(orientation_) {
	}

	species_id_t species_id;
	orientation_t orientation;
};


class reaction_t {
public:
	std::string name;
	float_t rate_constant;
	float_t min_noreaction_p;
	std::vector<species_with_orientation_t> reactants;
	std::vector<species_with_orientation_t> products;

	void dump(const std::string ind) const;
};

} // namespace mcell

#endif // SRC4_REACTION_H_
