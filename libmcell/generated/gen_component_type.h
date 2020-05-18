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

#ifndef API_GEN_COMPONENT_TYPE_H
#define API_GEN_COMPONENT_TYPE_H

#include "../api/common.h"
#include "../api/component_instance.h"

namespace MCell {
namespace API {

class ComponentInstance;

#define COMPONENT_TYPE_CTOR() \
    ComponentType( \
        const std::string& name_, \
        const std::vector<std::string> states_ = std::vector<std::string>() \
    ) { \
      class_name = "ComponentType"; \
      name = name_; \
      states = states_; \
      postprocess_in_ctor();\
      check_semantics();\
    }

class GenComponentType: public BaseDataClass {
public:
  void postprocess_in_ctor() override {}
  void check_semantics() const override;
  std::string to_str(const std::string ind="") const override;

  // --- attributes ---
  std::vector<std::string> states;
  virtual void set_states(const std::vector<std::string> new_states_) {
    states = new_states_;
  }
  virtual std::vector<std::string> get_states() const {
    return states;
  }

  // --- methods ---
  virtual ComponentInstance inst(const std::string& state = "STATE_UNSET", const int bond = BOND_UNBOUND) = 0;
  virtual ComponentInstance inst(const int state = STATE_UNSET_INT, const int bond = BOND_UNBOUND) = 0;
}; // GenComponentType

class ComponentType;
py::class_<ComponentType> define_pybinding_ComponentType(py::module& m);
} // namespace API
} // namespace MCell

#endif // API_GEN_COMPONENT_TYPE_H
