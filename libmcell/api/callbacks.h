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

#ifndef LIBMCELL_API_CALLBACKS_H_
#define LIBMCELL_API_CALLBACKS_H_

#include <functional>

#ifdef _WIN64
// fix for _hypot compilation issue
#define _hypot hypot
#include <cmath>
#endif
#include "pybind11/include/pybind11/pybind11.h"

#include "api/common.h"
#include "defines.h"
#include "bng/bng_defines.h"

namespace MCell {
namespace API {

class Model;
class MolWallHitInfo;
class ReactionInfo;

typedef std::function<void(std::shared_ptr<API::MolWallHitInfo>, pybind11::object)>
  wall_hit_callback_function_t;

typedef std::function<void(std::shared_ptr<API::ReactionInfo>, pybind11::object)>
  rxn_callback_function_t;

// not generated
// TODO: allow multiple callbacks
class Callbacks {
public:
  // model_ is nullptr in MDL mode
  Callbacks(Model* model_);

  Model* model;

  // -------------------- mol-wall hit callbacks --------------------
  void register_mol_wall_hit_callback(
      const wall_hit_callback_function_t func,
      py::object context,
      const geometry_object_id_t geometry_object_id,
      const species_id_t species_id
  ) {
    assert(model != nullptr);
    mol_wall_hit_callback_function = func;
    mol_wall_hit_context = context;
    mol_wall_hit_object_id = geometry_object_id;
    mol_wall_hit_species_id = species_id;
  }

  bool needs_callback_for_mol_wall_hit(
      const geometry_object_id_t geometry_object_id,
      const species_id_t species_id) {
    // model may be nullptr
    return mol_wall_hit_callback_function != nullptr &&
        (mol_wall_hit_object_id == GEOMETRY_OBJECT_ID_INVALID || mol_wall_hit_object_id == geometry_object_id) &&
        (mol_wall_hit_species_id == SPECIES_ID_INVALID || mol_wall_hit_species_id == species_id);
  }

  void do_mol_wall_hit_callback(std::shared_ptr<MolWallHitInfo> info);

  wall_hit_callback_function_t mol_wall_hit_callback_function;
  py::object mol_wall_hit_context;
  geometry_object_id_t mol_wall_hit_object_id; // GEOMETRY_OBJECT_ID_INVALID - any object may be hit
  species_id_t mol_wall_hit_species_id; // SPECIES_ID_INVALID - any species may be hit


  // -------------------- reaction callbacks --------------------
  void register_rxn_callback(
      const rxn_callback_function_t func,
      py::object context,
      const BNG::rxn_rule_id_t rxn_rule_id_
  ) {
    assert(model != nullptr);
    rxn_callback_function = func;
    rxn_context = context;
    rxn_rule_id = rxn_rule_id_;
  }

  bool needs_rxn_callback(
      const BNG::rxn_rule_id_t rxn_rule_id_) {
    return rxn_callback_function != nullptr &&
        rxn_rule_id == rxn_rule_id_;
  }

  void do_rxn_callback(std::shared_ptr<ReactionInfo> info);

  rxn_callback_function_t rxn_callback_function;
  py::object rxn_context;
  BNG::rxn_rule_id_t rxn_rule_id;

};

} /* namespace API */
} /* namespace MCell */

#endif /* LIBMCELL_API_CALLBACKS_H_ */
