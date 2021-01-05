/******************************************************************************
 *
 * Copyright (C) 2020 by
 * The Salk Institute for Biological Studies
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

#ifndef API_BASE_CHKPT_MOL_H
#define API_BASE_CHKPT_MOL_H

#include "generated/gen_base_chkpt_mol.h"
#include "api/common.h"

namespace MCell {

class Molecule;

namespace API {

// macro for GenChpt*Mol construction, used in Chpt*Mol constructors
#define GEN_CHKPT_MOL_CTOR(GEN_CLASS_NAME, M, ID_SPECIES_MAP, TIME_UNIT) \
    GEN_CLASS_NAME(\
        M.id, \
        ID_SPECIES_MAP.find(M.species_id)->second, \
        M.diffusion_time * TIME_UNIT, \
		    M.birthday * TIME_UNIT, \
		    ((M.unimol_rx_time != TIME_INVALID) ? M.unimol_rx_time * TIME_UNIT : FLT_UNSET) \
    )

class BaseChkptMol: public GenBaseChkptMol {
public:
  BASE_CHKPT_MOL_CTOR()


  bool export_as_string_without_newlines() const override {
    return true;
  }

  bool export_even_if_already_exported() const override {
    return true;
  }
};

} // namespace API
} // namespace MCell

#endif // API_BASE_CHKPT_MOL_H
