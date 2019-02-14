/*
 * partition.h
 *
 *  Created on: Feb 11, 2019
 *      Author: adam
 */

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
