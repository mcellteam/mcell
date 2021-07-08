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

#ifndef API_INITIAL_SURFACE_RELEASE_H
#define API_INITIAL_SURFACE_RELEASE_H

#include "generated/gen_initial_surface_release.h"
#include "api/api_common.h"
#include "api/complex.h"

namespace MCell {
namespace API {

class InitialSurfaceRelease: public GenInitialSurfaceRelease {
public:
  INITIAL_SURFACE_RELEASE_CTOR()

  void check_semantics() const override {
    if (get_num_set(number_to_release, density) != 1) {
      throw ValueError(S("Exactly one of ") + NAME_NUMBER_TO_RELEASE + " and " +
          NAME_DENSITY + " must be set.");
    }

    if (is_set(complex->compartment_name)) {
      throw ValueError(S("Compartment of the complex to be released in ") +
          NAME_CLASS_INITIAL_SURFACE_RELEASE + " must not be set.");
    }
  }

};

} // namespace API
} // namespace MCell

#endif // API_INITIAL_SURFACE_RELEASE_H
