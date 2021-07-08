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

#ifndef API_RELEASE_PATTERN_H
#define API_RELEASE_PATTERN_H

#include "generated/gen_release_pattern.h"
#include "api/api_common.h"

namespace MCell {
namespace API {

class ReleasePattern: public GenReleasePattern {
public:
  RELEASE_PATTERN_CTOR()

  void check_semantics() const override {
    GenReleasePattern::check_semantics();

    if (train_interval < train_duration) {
      throw ValueError(S(NAME_RELEASE_PATTERN) + name + ": " +
          NAME_TRAIN_INTERVAL + " is shorter than " + NAME_TRAIN_DURATION + ".");
    }
  }
};

} // namespace API
} // namespace MCell

#endif // API_RELEASE_PATTERN_H
