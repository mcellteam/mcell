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

// NOTE: this file should be called config.h, however with MSVC, there is
// an include collision and pybind11 includes it instead of some other file


#ifndef API_CONFIG_H
#define API_CONFIG_H

#include "generated/gen_config.h"
#include "api/api_common.h"

namespace MCell {
namespace API {

class Config: public GenConfig {
public:
  CONFIG_CTOR()

  void check_semantics() const override {
    GenConfig::check_semantics();
    if (cmp_gt(subpartition_dimension, partition_dimension, EPS)) {
      throw ValueError(S("Value ") + NAME_SUBPARTITION_DIMENSION + " must be smaller or equal than " + NAME_PARTITION_DIMENSION + ".");
    }

    if (is_set(initial_partition_origin)) {
      if (initial_partition_origin.size() != 3) {
        throw ValueError(S("Value ") + NAME_INITIAL_PARTITION_ORIGIN + " must be a vector of three floating point values.");
      }
    }
  }

};

} // namespace API
} // namespace MCell

#endif // API_CONFIG_H
