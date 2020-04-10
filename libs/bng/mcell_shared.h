/*
 * mcell_shared.h
 *
 *  Created on: Apr 1, 2020
 *      Author: ahusar
 */

#ifndef LIBS_BNG_MCELL_SHARED_H_
#define LIBS_BNG_MCELL_SHARED_H_

namespace MCell {
// same as in mcell_structs but renamed to make sure it is used correctly
enum species_flag_t {

  SPECIES_FLAG_ON_GRID = 0x01,
  SPECIES_FLAG_IS_SURFACE = 0x02,
  SPECIES_FLAG_CAN_VOLVOL = 0x10, // can vol vol react?
  SPECIES_FLAG_CAN_VOLSURF = 0x20,
  SPECIES_FLAG_CAN_SURFSURF = 0x80,
  SPECIES_FLAG_CANT_INITIATE = 0x400, // must not be set, not sure what to do with this yet (at least for some cases)
  SPECIES_FLAG_CAN_SURFSURFSURF = 0x20000, // 0x20000 - not supported
  SPECIES_FLAG_SET_MAX_STEP_LENGTH = 0x80000,
  SPECIES_FLAG_CAN_REGION_BORDER = 0x100000, // CAN_REGION_BORDER, // 0x100000
  SPECIES_FLAG_EXTERNAL_SPECIES = 0x400000 // 0x400000 - not supported
};

}

#endif /* LIBS_BNG_MCELL_SHARED_H_ */
