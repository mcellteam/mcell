/******************************************************************************
 *
 * Copyright (C) 2021 by
 * The Salk Institute for Biological Studies and
 * Pittsburgh Supercomputing Center, Carnegie Mellon University
 *
 * Use of this source code is governed by an MIT-style
 * license that can be found in the LICENSE file or at
 * https://opensource.org/licenses/MIT.
 *
******************************************************************************/

#ifndef SRC4_REGION_UTIL_H_
#define SRC4_REGION_UTIL_H_

#include "defines.h"

namespace MCell {

class Partition;
class Molecule;

namespace RegionUtils {

/* ALL_INSIDE: flag that indicates that all reactants lie inside their
 *             respective restrictive regions
 * ALL_OUTSIDE: flag that indicates that all reactants lie outside
 *              their respective restrictive regions
 * SURF1_IN_SURF2_OUT: flag that indicates that  reactant "sm_1" lies
 *                     inside and reactant "sm_2" lies outside of their
 *                     respective restrictive regions
 * SURF1_OUT_SURF2_IN: flag that indicates that  reactant "sm_1" lies outside
 *                     and reactant "sm_2" lies inside of their
 *                     respective restrictive regions
 * SURF1_IN: flag that indicates that only reactant "sm_1" has
 *          restrictive regions on the object and it lies
 *          inside its restrictive region.
 * SURF1_OUT: flag that indicates that only reactant "sm_1" has
 *            restrictive regions on the object and it lies
 *            outside its restrictive region.
 * SURF2_IN: flag that indicates that only reactant "sm_2" has
 *           restrictive regions on the object and it lies
 *           inside its restrictive region.
 * SURF2_OUT: flag that indicates that only reactant "sm_2" has
 *            restrictive regions on the object and it lies
 *            outside its restrictive region.  */
enum surf_mol_region_relationship_flag_t {
  ALL_INSIDE = 0x01,
  ALL_OUTSIDE = 0x02,
  SURF1_IN_SURF2_OUT = 0x04,
  SURF1_OUT_SURF2_IN = 0x08,
  SURF1_IN = 0x10,
  SURF1_OUT = 0x20,
  SURF2_IN = 0x40,
  SURF2_OUT = 0x80
};


// returns flags composed of surf_mol_region_relationship_flag_t
uint determine_molecule_region_topology(
    Partition& p,
    const Molecule* reacA,
    const Molecule* reacB,
    const bool is_unimol,
    RegionIndicesSet& rlp_wall_1,
    RegionIndicesSet& rlp_wall_2,
    RegionIndicesSet& rlp_obj_1,
    RegionIndicesSet& rlp_obj_2);


bool product_tile_can_be_reached(
    const Partition& p,
    const wall_index_t wall_index,
    const bool is_unimol,
    const uint sm_bitmask,
    const RegionIndicesSet& rlp_wall_1,
    const RegionIndicesSet& rlp_wall_2,
    const RegionIndicesSet& rlp_obj_1,
    const RegionIndicesSet& rlp_obj_2);

} // namespace RegionUtil
} // namespace MCell

#endif // SRC4_REGION_UTIL_H_
