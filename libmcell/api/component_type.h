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

#ifndef API_COMPONENT_TYPE_H
#define API_COMPONENT_TYPE_H

#include "generated/gen_component_type.h"
#include "api/api_common.h"
#include "api/component.h"

namespace MCell {
namespace API {

class ComponentType: public GenComponentType, public std::enable_shared_from_this<ComponentType> {
public:
  COMPONENT_TYPE_CTOR()

  std::shared_ptr<Component> inst(const std::string& state, const int bond) override {
    return std::make_shared<Component>(shared_from_this(), state, bond);
  }

  std::shared_ptr<Component> inst(const int state, const int bond) override {
    return std::make_shared<Component>(shared_from_this(), std::to_string(state), bond);
  }

  bool __eq__(const ComponentType& other) const override;

  std::string to_bngl_str() const override;

  // added manually
  std::string get_canonical_name() const;

  // needed when defining a set of ComponentTypes
  bool operator < (const ComponentType& other) const;
};

} // namespace API
} // namespace MCell

#endif // API_COMPONENT_TYPE_H
