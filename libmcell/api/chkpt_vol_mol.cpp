/******************************************************************************
 *
 * Copyright (C) 2020 by
 * The Salk Institute for Biological Studies
 *
 * Use of this source code is governed by an MIT-style
 * license that can be found in the LICENSE file or at
 * https://opensource.org/licenses/MIT.
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
    const IdSpeciesMap& id_species_map, const double time_unit,
    const double length_unit) :
    GEN_CHKPT_MOL_CTOR(GenChkptVolMol, vm, id_species_map, time_unit) {

  assert(vm.id != MOLECULE_ID_INVALID);

  type = MoleculeType::VOLUME;
  pos = vm.v.pos * Vec3(length_unit);
}


} // namespace API
} // namespace MCell
