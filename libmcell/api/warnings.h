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

#ifndef API_WARNINGS_H
#define API_WARNINGS_H

#include "generated/gen_warnings.h"
#include "api/api_common.h"

namespace MCell {
namespace API {

class Warnings: public GenWarnings {
public:
  WARNINGS_CTOR()

  // must be called also manually during model initialization
  void check_semantics() const override {
    if (high_reaction_probability == WarningLevel::ERROR) {
      throw ValueError(S(NAME_CLASS_WARNINGS) + "." + NAME_HIGH_REACTION_PROBABILITY + " must not be set to " +
          NAME_ENUM_WARNING_LEVEL + "." + NAME_EV_ERROR + ".");
    }
  }
};

} // namespace API
} // namespace MCell

#endif // API_WARNINGS_H
