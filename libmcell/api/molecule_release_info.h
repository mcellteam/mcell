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

#ifndef API_MOLECULE_RELEASE_INFO_H
#define API_MOLECULE_RELEASE_INFO_H

#include "generated/gen_molecule_release_info.h"
#include "api/api_common.h"

namespace MCell {
namespace API {

class MoleculeReleaseInfo: public GenMoleculeReleaseInfo {
public:
  MOLECULE_RELEASE_INFO_CTOR()

  void check_semantics() const override {
    if (location.size() != 3) {
      throw ValueError(S("Value of ") + NAME_LOCATION + " must be a triplet of floats.");
    }
  }
};

} // namespace API
} // namespace MCell

#endif // API_MOLECULE_RELEASE_INFO_H
