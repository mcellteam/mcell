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

#ifndef API_MOL_WALL_HIT_INFO_H
#define API_MOL_WALL_HIT_INFO_H

#include "generated/gen_mol_wall_hit_info.h"
#include "api/api_common.h"

namespace MCell {
namespace API {

class MolWallHitInfo: public GenMolWallHitInfo {
public:
  MolWallHitInfo() {
  }
  MolWallHitInfo(DefaultCtorArgType) {
  }

  void dump() {
    std::cout << to_str();
  }

  // extra information to be converted in Callbacks
  geometry_object_id_t geometry_object_id; // to geometry_object
  wall_index_t partition_wall_index; // to wall_index
};

} // namespace API
} // namespace MCell

#endif // API_MOL_WALL_HIT_INFO_H
