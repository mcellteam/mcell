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

#ifndef LIBMCELL_API_BASE_CHKPT_MOL_H
#define LIBMCELL_API_BASE_CHKPT_MOL_H

#include "generated/gen_base_chkpt_mol.h"
#include "api/api_common.h"

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
		    M.flags, \
		    ((M.unimol_rx_time != TIME_INVALID) ? M.unimol_rx_time * TIME_UNIT : FLT_UNSET) \
    )

class BaseChkptMol: public GenBaseChkptMol {
public:
  BASE_CHKPT_MOL_CTOR()
  BaseChkptMol() {
  }

  void postprocess_in_ctor() override {
    set_all_custom_attributes_to_default();
  }

  void set_all_custom_attributes_to_default() override {
    type = MoleculeType::UNSET;
  }

  bool export_as_string_without_newlines() const override {
    return true;
  }

  bool export_even_if_already_exported() const override {
    return true;
  }

  // used when casting to a derived class
  MoleculeType type;
};

} // namespace API
} // namespace MCell

#endif // LIBMCELL_API_BASE_CHKPT_MOL_H
