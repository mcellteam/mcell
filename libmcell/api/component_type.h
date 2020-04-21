/******************************************************************************
 *
 * Copyright (C) 2020 by
 * The Salk Institute for Biological Studies and
 * Pittsburgh Supercomputing Center, Carnegie Mellon University
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

#include <string>

#include "../generated/gen_component_type.h"
#include "common.h"
#include "component_instance.h"

namespace MCell {
namespace API {

class ComponentType: public GenComponentType {
public:
  COMPONENT_TYPE_CTOR()

  ComponentInstance inst(const std::string& state, const int bond) override {
    return ComponentInstance(this, state, bond);
  }

  ComponentInstance inst(const int state, const int bond) override {
    return ComponentInstance(this, std::to_string(state), bond);
  }
};


} // namespace API
} // namespace MCell

#endif // API_COMPONENT_TYPE_H
