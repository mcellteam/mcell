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

#ifndef API_CHKPT_VOL_MOL_H
#define API_CHKPT_VOL_MOL_H

#include "generated/gen_chkpt_vol_mol.h"
#include "api/api_common.h"
#include "api/base_chkpt_mol.h"


namespace MCell {

class Molecule;

namespace API {

class ChkptVolMol: public GenChkptVolMol {
public:
  CHKPT_VOL_MOL_CTOR()

  ChkptVolMol(
    const MCell::Molecule& vm,
    const IdSpeciesMap& id_species_map, const double time_unit,
    const double length_unit);

  void postprocess_in_ctor() override {
    set_all_custom_attributes_to_default();
  }

  void set_all_custom_attributes_to_default() override {
    type = MoleculeType::VOLUME;
  }

  virtual ~ChkptVolMol() {
  };
};

} // namespace API
} // namespace MCell

#endif // API_CHKPT_VOL_MOL_H
