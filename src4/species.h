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

#ifndef SRC4_SPECIES_H_
#define SRC4_SPECIES_H_

#include <string>
#include "defines.h"

namespace mcell {

class species_t {
public:
	species_id_t species_id;

	uint32_t mcell3_species_id;
	float_t D; // diffusion constant
	std::string name;
	float_t space_step;
	float_t time_step; // in standard time

	void dump(const std::string ind);
};

} /* namespace mcell */

#endif /* SRC4_SPECIES_H_ */
