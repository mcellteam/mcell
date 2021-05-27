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
