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

#ifndef API_COMPONENT_INSTANCE_H
#define API_COMPONENT_INSTANCE_H

#include "generated/gen_component.h"
#include "api/api_common.h"

namespace MCell {
namespace API {

class Component: public GenComponent {
public:
  COMPONENT_CTOR()

  // default __eq__ operator is sufficient

  // needed when defining a set of ComponentInstances
  bool operator < (const Component& other) const;

  std::string to_bngl_str() const override;
};

} // namespace API
} // namespace MCell

#endif // API_COMPONENT_INSTANCE_H
