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

// NOTE: this file should be called config.h, however with MSVC, there is
// an include collision and pybind11 includes it instead of some other file


#ifndef LIBMCELL_API_CONFIG_H
#define LIBMCELL_API_CONFIG_H

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

#endif // LIBMCELL_API_CONFIG_H
