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

#ifndef API_CHKPT_SURF_MOL_H
#define API_CHKPT_SURF_MOL_H

#include "generated/gen_chkpt_surf_mol.h"
#include "api/api_common.h"
#include "api/base_chkpt_mol.h"


namespace MCell {

class Partition;
class Molecule;

namespace API {

class ChkptSurfMol: public GenChkptSurfMol {
public:
  CHKPT_SURF_MOL_CTOR()

  ChkptSurfMol(
      const MCell::Molecule& sm,
      const IdSpeciesMap& id_species_map, const double time_unit,
      const double length_unit,
      const MCell::Partition& p,
      const IdGeometryObjectMap& id_geometry_object_map);

  void postprocess_in_ctor() override {
    set_all_custom_attributes_to_default();
  }

  void set_all_custom_attributes_to_default() override {
    type = MoleculeType::SURFACE;
  }

  virtual ~ChkptSurfMol() {
  }
};

} // namespace API
} // namespace MCell

#endif // API_CHKPT_SURF_MOL_H
