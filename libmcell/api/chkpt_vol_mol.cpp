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

#include "api/chkpt_vol_mol.h"
#include "api/base_chkpt_mol.h"

#include "src4/molecule.h"

using namespace std;

namespace MCell {
namespace API {

ChkptVolMol::ChkptVolMol(
    const MCell::Molecule& vm,
    const IdSpeciesMap& id_species_map, const float_t time_unit,
    const float_t length_unit) :
    GEN_CHKPT_MOL_CTOR(GenChkptVolMol, vm, id_species_map, time_unit) {

  assert(vm.id != MOLECULE_ID_INVALID);

  pos = vm.v.pos * Vec3(length_unit);
}


} // namespace API
} // namespace MCell
