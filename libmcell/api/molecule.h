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

#ifndef API_MOLECULE_H
#define API_MOLECULE_H

#include "generated/gen_molecule.h"
#include "api/api_common.h"
#include "api/species.h"

namespace MCell {
namespace API {

class Molecule: public GenMolecule {
public:
  MOLECULE_CTOR_NOARGS()

  void remove() override;
};

} // namespace API
} // namespace MCell

#endif // API_MOLECULE_H
