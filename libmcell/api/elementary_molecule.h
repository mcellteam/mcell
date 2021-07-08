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

#ifndef API_ELEMENTARY_MOLECULE_INSTANCE_H
#define API_ELEMENTARY_MOLECULE_INSTANCE_H

#include "generated/gen_elementary_molecule.h"
#include "api/api_common.h"

namespace MCell {
namespace API {

class ElementaryMolecule: public GenElementaryMolecule {
public:
  ELEMENTARY_MOLECULE_CTOR()

  bool __eq__(const ElementaryMolecule& other) const override;

  std::string to_bngl_str(const bool with_compartment = true) const override;

  bool is_surf() const;
};

} // namespace API
} // namespace MCell

#endif // API_ELEMENTARY_MOLECULE_INSTANCE_H
