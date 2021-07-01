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

#include "api/chkpt_surf_mol.h"

#include "api/api_utils.h"
#include "api/geometry_object.h"
#include "src4/molecule.h"
#include "src4/geometry.h"
#include "src4/partition.h"

using namespace std;

namespace MCell {
namespace API {

ChkptSurfMol::ChkptSurfMol(
    const MCell::Molecule& sm,
    const IdSpeciesMap& id_species_map, const double time_unit,
    const double length_unit,
    const MCell::Partition& p,
    const IdGeometryObjectMap& id_geometry_object_map) :
    GEN_CHKPT_MOL_CTOR(GenChkptSurfMol, sm, id_species_map, time_unit) {

  assert(sm.id != MOLECULE_ID_INVALID);

  type = MoleculeType::SURFACE;

  pos = sm.s.pos * Vec2(length_unit);
  orientation = convert_mcell_orientation(sm.s.orientation);

  const Wall& w = p.get_wall(sm.s.wall_index);

  auto it = id_geometry_object_map.find(w.object_id);
  assert(it != id_geometry_object_map.end());
  geometry_object = it->second;

  wall_index = geometry_object->get_object_wall_index(sm.s.wall_index);
  grid_tile_index = sm.s.grid_tile_index;
}

} // namespace API
} // namespace MCell
