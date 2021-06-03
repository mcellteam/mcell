/******************************************************************************
 *
 * Copyright (C) 2020-2021 by
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
#include "molecule.h"

using namespace std;

namespace MCell {
namespace API {

Callbacks::Callbacks(Model* model_)
  : model(model_) {
  assert(model != nullptr);
}


void Callbacks::register_mol_wall_hit_callback(
    const mol_wall_hit_callback_function_t func,
    py::object context,
    const geometry_object_id_t geometry_object_id,
    const BNG::species_id_t species_id
) {
  assert(model != nullptr);

  auto it_geom_obj = mol_wall_hit_callbacks.find(geometry_object_id);
  if (it_geom_obj != mol_wall_hit_callbacks.end() && it_geom_obj->second.count(species_id) != 0) {

    string geom_name;
    if (geometry_object_id != GEOMETRY_OBJECT_ID_INVALID) {
      geom_name = model->get_world()->get_geometry_object(geometry_object_id).name;
    }
    else {
      geom_name = "any";
    }

    string species_name;
    if (species_id != BNG::SPECIES_ID_INVALID) {
      species_name = model->get_world()->get_all_species().get(species_id).name;
    }
    else {
      species_name = "any";
    }

    throw RuntimeError(S("Cannot register two callbacks for an identical pair or geometry object and species id, ") +
        " error while trying to register second callback for geometry object '" + geom_name + "' and species '" + species_name +"'.");
  }

  mol_wall_hit_callbacks[geometry_object_id][species_id] = MolWallHitCallbackInfo(func, context, geometry_object_id, species_id);

  // make sure that the species_id won't change in the future
  if (species_id != BNG::SPECIES_ID_INVALID) {
    model->get_world()->get_all_species().get(species_id).clear_flag(BNG::SPECIES_FLAG_IS_REMOVABLE);
  }
}


void Callbacks::do_mol_wall_hit_callbacks(std::shared_ptr<MolWallHitInfo> info) {

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

  geometry_object_id_t geometry_object_id = info->geometry_object->geometry_object_id;

  // only one partition for now
  const MCell::Molecule& m = model->get_world()->get_partition(PARTITION_ID_INITIAL).get_m(info->molecule_id);

  // call callback for all matching registered callbacks
  auto it_specific_geom_obj = mol_wall_hit_callbacks.find(geometry_object_id);
  if (it_specific_geom_obj != mol_wall_hit_callbacks.end()) {
    do_mol_wall_hit_callback_for_specific_and_any_species(
        info, m.species_id, it_specific_geom_obj->second);
  }

  auto it_any_geom_obj = mol_wall_hit_callbacks.find(GEOMETRY_OBJECT_ID_INVALID);
  if (it_any_geom_obj != mol_wall_hit_callbacks.end()) {
    do_mol_wall_hit_callback_for_specific_and_any_species(
        info, m.species_id, it_any_geom_obj->second);
  }
}


void Callbacks::do_mol_wall_hit_callback_for_specific_and_any_species(
    std::shared_ptr<MolWallHitInfo> info,
    const BNG::species_id_t specific_species_id,
    const SpeciesMolWallHitCallbackInfoMap& species_map) {

  auto it_specific_species = species_map.find(specific_species_id);
  if (it_specific_species != species_map.end()) {
    do_individual_mol_wall_hit_callback(info, it_specific_species->second);
  }

  auto it_any_species = species_map.find(BNG::SPECIES_ID_INVALID);
  if (it_any_species != species_map.end()) {
    do_individual_mol_wall_hit_callback(info, it_any_species->second);
  }
}


void Callbacks::do_individual_mol_wall_hit_callback(
    std::shared_ptr<MolWallHitInfo> info,
    MolWallHitCallbackInfo callback_function_and_context) {

  // acquire GIL before calling Python code
  py::gil_scoped_acquire acquire;

  // call the actual callback
  callback_function_and_context.callback_function(info, callback_function_and_context.context);
}


void Callbacks::register_rxn_callback(
    const rxn_callback_function_t func,
    py::object context,
    const BNG::rxn_rule_id_t rxn_rule_id
) {
  assert(model != nullptr);
  assert(rxn_rule_id != BNG::RXN_RULE_ID_INVALID);

  if (rxn_callbacks.count(rxn_rule_id) != 0) {
    std::string name = model->get_world()->get_all_rxns().get(rxn_rule_id)->to_str();
    throw RuntimeError(S("Each reaction rule can have only a single callback, error while trying to register ") +
        "second callback for " + name + ".");
  }

  rxn_callbacks[rxn_rule_id] = RxnCallbackInfo(func, context, rxn_rule_id);
}


void Callbacks::do_rxn_callback(std::shared_ptr<ReactionInfo> info) {
  // select the correct callback
  assert(rxn_callbacks.count(info->rxn_rule_id) != 0);
  const RxnCallbackInfo& specific_callback = rxn_callbacks[info->rxn_rule_id];

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

  // acquire GIL before calling Python code
  py::gil_scoped_acquire acquire;

  // call the actual callback
  specific_callback.callback_function(info, specific_callback.context);
}

} /* namespace API */
} /* namespace MCell */
