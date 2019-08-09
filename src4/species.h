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

#include "mcell_structs.h"

namespace mcell {

// same as in mcell_structs but renamed to make sure it is used correctly
enum species_flag_t {

  SPECIES_FLAG_ON_GRID = ON_GRID, // 0x01
  SPECIES_FLAG_CAN_VOLVOL = CAN_VOLVOL, // 0x10  - can vol vol react?
  SPECIES_FLAG_CAN_VOLSURF = CAN_VOLSURF, // 0x20
  SPECIES_FLAG_CAN_SURFSURF = CAN_SURFSURF, // 0x80
  SPECIES_FLAG_CANT_INITIATE = CANT_INITIATE, // 0x400 - not sure what to do with this yet
  SPECIES_FLAG_CAN_SURFSURFSURF = CAN_SURFSURFSURF, // 0x20000 - not supported
  SPECIES_FLAG_SET_MAX_STEP_LENGTH = SET_MAX_STEP_LENGTH, // 0x80000
  SPECIES_FLAG_CAN_REGION_BORDER = CAN_REGION_BORDER, // 0x100000 - not supported
  SPECIES_FLAG_EXTERNAL_SPECIES = EXTERNAL_SPECIES // 0x400000 - not supported
};
/**
 * Holds information on one species type.
 */
class species_t {
public:
  species_id_t species_id;

  uint mcell3_species_id;
  float_t D; // diffusion constant
  std::string name;
  float_t space_step;
  float_t time_step; // in standard time

  uint flags;

  bool has_flag(species_flag_t flag) const {
    return (flags & flag) != 0;
  }

  bool is_surf() const {
    return has_flag(SPECIES_FLAG_ON_GRID);
  }

  bool is_vol() const {
    return !has_flag(SPECIES_FLAG_ON_GRID);
  }

  bool can_diffuse() const {
    return D != DIFFUSION_CONSTANT_ZER0;
  }

  void dump(const std::string ind) const;
  static void dump_array(const std::vector<species_t>& vec);
};

} // namespace mcell

#endif // SRC4_SPECIES_H_

