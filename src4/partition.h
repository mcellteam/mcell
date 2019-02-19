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

#ifndef SRC4_PARTITION_H_
#define SRC4_PARTITION_H_

#include "defines.h"
#include "molecule.h"
#include "molecule.h"

namespace mcell {

class partition_t {
	// left, bottom, closest (lowest z) point of the partition
  vec3_t origin;

  // vector containing all volume molecules in this partition
  std::vector< /* molecule index*/ volume_molecule_t> volume_molecules;

  // arrays of indices to the volume_molecules array where each array corresponds to a given time step
  std::vector< /* diffusion time step index */
		std::pair< float_t, std::vector< molecule_index_t > > > volume_molecule_indices_per_time_step;

  // TBD
  //std::vector< /* subpartition index */
  //std::vector < /* diffusion time step index */ subpartition_mask_t > > volume_molecules_subpartition_masks;

  //TBD: std::vector< /* surface molecule index */ surface_molecule> surface_molecules;
  //TBD: std::vector< /* subpartition index */ subpartition_mask > surface_molecules_subpatition_masks;

};

} // namespace mcell

#endif /* SRC4_PARTITION_H_ */
