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

#include "api/callbacks.h"

#include "api/model.h"
#include "api/geometry_object.h"

namespace MCell {
namespace API {

Callbacks::Callbacks(Model* model_)
  : model(model_),
    mol_wall_hit_callback_function(nullptr),
    mol_wall_hit_object_id(GEOMETRY_OBJECT_ID_INVALID),
    mol_wall_hit_species_id(SPECIES_ID_INVALID)
  {
  assert(model != nullptr);
}


void Callbacks::do_mol_wall_hit_callback(std::shared_ptr<MolWallHitInfo> info) {
  info->geometry_object = model->get_geometry_object_with_id(info->geometry_object_id);
  assert(is_set(info->geometry_object));
  assert(info->partition_wall_index >= info->geometry_object->first_wall_index);
  info->wall_index = info->partition_wall_index - info->geometry_object->first_wall_index;
  // call the actual callback
  mol_wall_hit_callback_function(info, mol_wall_hit_context);
}

} /* namespace API */
} /* namespace MCell */
