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
