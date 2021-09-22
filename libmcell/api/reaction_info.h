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

#ifndef API_VOLUME_REACTION_INFO_H
#define API_VOLUME_REACTION_INFO_H

#include "generated/gen_reaction_info.h"
#include "api/api_common.h"
#include "bng/bng_defines.h"

namespace MCell {
namespace API {

class ReactionInfo: public GenReactionInfo {
public:
  void dump() {
    std::cout << to_str();
  }

  ReactionInfo() {
    set_all_custom_attributes_to_default();
  }

  ReactionInfo(DefaultCtorArgType) {
    set_all_custom_attributes_to_default();
  }

  void set_all_custom_attributes_to_default() {
    // setting all (also inherited) attributes

    rxn_rule_id = BNG::RXN_RULE_ID_INVALID;
    geometry_object_id = GEOMETRY_OBJECT_ID_INVALID;
    partition_wall_index = WALL_INDEX_INVALID;

    // inherited attributes
    type = ReactionType::UNSET;
    reaction_rule = nullptr;
    time = FLT_UNSET;
    pos3d.clear();
    geometry_object = nullptr;
    wall_index = -1;
    pos2d.clear();
  }

  // extra information to be converted in Callbacks
  BNG::rxn_rule_id_t rxn_rule_id;

  geometry_object_id_t geometry_object_id; // to geometry_object
  wall_index_t partition_wall_index;
};

} // namespace API
} // namespace MCell

#endif // API_VOLUME_REACTION_INFO_H
