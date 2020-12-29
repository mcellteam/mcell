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

#ifndef API_CHKPT_SURF_MOL_H
#define API_CHKPT_SURF_MOL_H

#include "generated/gen_chkpt_surf_mol.h"
#include "api/common.h"
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
      const IdSpeciesMap& id_species_map, const float_t time_unit,
      const float_t length_unit,
      const MCell::Partition& p,
      const IdGeometryObjectMap& id_geometry_object_map
  );

  virtual ~ChkptSurfMol() {
  }
};

} // namespace API
} // namespace MCell

#endif // API_CHKPT_SURF_MOL_H
