/******************************************************************************
 *
 * Copyright (C) 2020-2021 by
 * The Salk Institute for Biological Studies
 *
 * Use of this source code is governed by an MIT-style
 * license that can be found in the LICENSE file or at
 * https://opensource.org/licenses/MIT.
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
#ifdef _MSC_VER
#undef HAVE_UNISTD_H
#undef HAVE_SYS_TIME_H
#endif
#include "pybind11/include/pybind11/pybind11.h"

#include "api/api_common.h"
#include "defines.h"
#include "bng/bng_defines.h"

namespace MCell {
namespace API {

class Model;
class MolWallHitInfo;
class ReactionInfo;

typedef std::function<void(std::shared_ptr<API::MolWallHitInfo>, pybind11::object)>
  mol_wall_hit_callback_function_t;

typedef std::function<bool(std::shared_ptr<API::ReactionInfo>, pybind11::object)>
  rxn_callback_function_t;


struct RxnCallbackInfo {
  RxnCallbackInfo() :
    callback_function(nullptr),
    rxn_rule_id(BNG::RXN_RULE_ID_INVALID) {
  }

  RxnCallbackInfo(
      const rxn_callback_function_t callback_function_,
      const py::object context_,
      const BNG::rxn_rule_id_t rxn_rule_id_) :
    callback_function(callback_function_),
    context(context_),
    rxn_rule_id(rxn_rule_id_) {
  }

  rxn_callback_function_t callback_function;
  py::object context;
  BNG::rxn_rule_id_t rxn_rule_id;
};


struct MolWallHitCallbackInfo {
  MolWallHitCallbackInfo() :
    callback_function(nullptr),
    geometry_object_id(GEOMETRY_OBJECT_ID_INVALID),
    species_id(BNG::SPECIES_ID_INVALID) {
  }

  MolWallHitCallbackInfo(
      const mol_wall_hit_callback_function_t callback_function_,
      const py::object context_,
      const geometry_object_id_t geometry_object_id_,
      const BNG::species_id_t species_id_) :
    callback_function(callback_function_),
    context(context_),
    geometry_object_id(geometry_object_id_),
    species_id(species_id_) {
  }

  mol_wall_hit_callback_function_t callback_function;
  py::object context;
  geometry_object_id_t geometry_object_id; // GEOMETRY_OBJECT_ID_INVALID - any object may be hit
  species_id_t species_id; // SPECIES_ID_INVALID - any species may be hit
};


// not generated
class Callbacks {
public:
  // model_ is nullptr in MDL mode
  Callbacks(Model* model_);

  Model* model;

  // -------------------- mol-wall hit callbacks --------------------
  void register_mol_wall_hit_callback(
      const mol_wall_hit_callback_function_t func,
      py::object context,
      const geometry_object_id_t geometry_object_id,
      const species_id_t species_id);

  bool needs_callback_for_mol_wall_hit(
      const geometry_object_id_t geometry_object_id,
      const species_id_t species_id) const {

    if (mol_wall_hit_callbacks.empty()) {
      return false;
    }

    auto it_specific_geom_obj = mol_wall_hit_callbacks.find(geometry_object_id);
    if (it_specific_geom_obj != mol_wall_hit_callbacks.end() &&
        needs_callback_for_mol_wall_hit(it_specific_geom_obj->second, species_id)) {
      return true;
    }

    auto it_any_geom_obj = mol_wall_hit_callbacks.find(GEOMETRY_OBJECT_ID_INVALID);
    if (it_any_geom_obj != mol_wall_hit_callbacks.end() &&
        needs_callback_for_mol_wall_hit(it_any_geom_obj->second, species_id)) {
      return true;
    }

    // no match
    return false;
  }

  void do_mol_wall_hit_callbacks(std::shared_ptr<MolWallHitInfo> info);

  // -------------------- reaction callbacks --------------------
  void register_rxn_callback(
      const rxn_callback_function_t func,
      py::object context,
      const BNG::rxn_rule_id_t rxn_rule_id
  );

  bool needs_rxn_callback(
      const BNG::rxn_rule_id_t rxn_rule_id) const {
    if (rxn_callbacks.empty()) {
      return false;
    }
    else {
      return rxn_callbacks.count(rxn_rule_id) != 0;
    }
  }

  // returns true if reaction should be cancelled
  bool do_rxn_callback(std::shared_ptr<ReactionInfo> info);

private:

  std::map<BNG::rxn_rule_id_t, RxnCallbackInfo> rxn_callbacks;

  typedef std::map<BNG::species_id_t, MolWallHitCallbackInfo> SpeciesMolWallHitCallbackInfoMap;
  std::map<geometry_object_id_t, SpeciesMolWallHitCallbackInfoMap> mol_wall_hit_callbacks;

  bool needs_callback_for_mol_wall_hit(
      const SpeciesMolWallHitCallbackInfoMap& species_info_map,
      const BNG::species_id_t species_id) const {

    if (species_info_map.count(species_id) != 0) {
      return true;
    }
    if (species_info_map.count(BNG::SPECIES_ID_INVALID) != 0) {
      return true;
    }
    return false;
  }

  void do_mol_wall_hit_callback_for_specific_and_any_species(
      std::shared_ptr<MolWallHitInfo> info,
      const BNG::species_id_t specific_species_id,
      const SpeciesMolWallHitCallbackInfoMap& species_map);

  void do_individual_mol_wall_hit_callback(
      std::shared_ptr<MolWallHitInfo> info, MolWallHitCallbackInfo callback_function_and_context);
};

} /* namespace API */
} /* namespace MCell */

#endif /* LIBMCELL_API_CALLBACKS_H_ */
