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

#ifndef API_COMPONENT_TYPE_H
#define API_COMPONENT_TYPE_H

#include "generated/gen_component_type.h"
#include "api/common.h"
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
