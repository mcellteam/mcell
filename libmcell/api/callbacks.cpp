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

#include "world.h"

namespace MCell {
namespace API {

Callbacks::Callbacks(Model* model_)
  : model(model_),

    mol_wall_hit_callback_function(nullptr),
    mol_wall_hit_object_id(GEOMETRY_OBJECT_ID_INVALID),
    mol_wall_hit_species_id(SPECIES_ID_INVALID),

    rxn_callback_function(nullptr),
    rxn_rule_id(BNG::RXN_RULE_ID_INVALID)
  {
  assert(model != nullptr);
}


void Callbacks::do_mol_wall_hit_callback(std::shared_ptr<MolWallHitInfo> info) {

  // set geometry data
  info->geometry_object = model->get_geometry_object_with_id(info->geometry_object_id);
  assert(is_set(info->geometry_object));
  assert(info->partition_wall_index >= info->geometry_object->first_wall_index);
  info->wall_index = info->partition_wall_index - info->geometry_object->first_wall_index;

  // convert units
  assert(model->get_world() != nullptr);
  info->time = info->time * model->get_world()->config.time_unit;
  info->pos3d = info->pos3d * Vec3(model->get_world()->config.length_unit);
  info->time_before_hit = info->time_before_hit * model->get_world()->config.time_unit;
  info->pos3d_before_hit = info->pos3d_before_hit * Vec3(model->get_world()->config.length_unit);

  // call the actual callback
  mol_wall_hit_callback_function(info, mol_wall_hit_context);
}


void Callbacks::do_rxn_callback(std::shared_ptr<ReactionInfo> info) {
  // set reaction rule object
  info->reaction_rule = model->get_reaction_rule_with_fwd_id(info->rxn_rule_id);
  assert(is_set(info->reaction_rule));

  // convert units
  assert(model->get_world() != nullptr);
  info->time = info->time * model->get_world()->config.time_unit;
  info->pos3d = info->pos3d * Vec3(model->get_world()->config.length_unit);

  const BNG::RxnRule* rxn = model->get_world()->get_all_rxns().get(info->rxn_rule_id);

  if (info->geometry_object_id != GEOMETRY_OBJECT_ID_INVALID) {
    info->geometry_object = model->get_geometry_object_with_id(info->geometry_object_id);
    assert(is_set(info->geometry_object));
    assert(info->partition_wall_index >= info->geometry_object->first_wall_index);
    info->wall_index = info->geometry_object->get_object_wall_index(info->partition_wall_index);
    info->pos2d = info->pos2d * Vec2(model->get_world()->config.length_unit);
  }

  if (info->geometry_object_id_surf_reac2 != GEOMETRY_OBJECT_ID_INVALID) {
    info->geometry_object_surf_reac2 = model->get_geometry_object_with_id(info->geometry_object_id_surf_reac2);
    assert(is_set(info->geometry_object_surf_reac2));
    assert(info->partition_wall_index_surf_reac2 >= info->geometry_object_surf_reac2->first_wall_index);
    info->wall_index_surf_reac2 = info->geometry_object_surf_reac2->get_object_wall_index(info->partition_wall_index_surf_reac2);
    info->pos2d_surf_reac2 = info->pos2d_surf_reac2 * Vec2(model->get_world()->config.length_unit);
  }

  // call the actual callback
  rxn_callback_function(info, rxn_context);
}

} /* namespace API */
} /* namespace MCell */
