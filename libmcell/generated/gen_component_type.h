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

#include "api/api_common.h"
#include "api/base_data_class.h"

namespace MCell {
namespace API {

class ComponentType;
class Component;
class PythonExportContext;

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
  void set_initialized() override;
  void set_all_attributes_as_default_or_unset() override;

  virtual bool __eq__(const ComponentType& other) const;
  virtual bool eq_nonarray_attributes(const ComponentType& other, const bool ignore_name = false) const;
  bool operator == (const ComponentType& other) const { return __eq__(other);}
  bool operator != (const ComponentType& other) const { return !__eq__(other);}
  std::string to_str(const std::string ind="") const override;

  std::string export_to_python(std::ostream& out, PythonExportContext& ctx) override;
  virtual std::string export_vec_states(std::ostream& out, PythonExportContext& ctx, const std::string& parent_name);


  // --- attributes ---
  std::vector<std::string> states;
  virtual void set_states(const std::vector<std::string> new_states_) {
    if (initialized) {
      throw RuntimeError("Value 'states' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    states = new_states_;
  }
  virtual std::vector<std::string>& get_states() {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return states;
  }

  // --- methods ---
  virtual std::shared_ptr<Component> inst(const std::string& state = "STATE_UNSET", const int bond = BOND_UNBOUND) = 0;
  virtual std::shared_ptr<Component> inst(const int state = STATE_UNSET_INT, const int bond = BOND_UNBOUND) = 0;
  virtual std::string to_bngl_str() const = 0;
}; // GenComponentType

class ComponentType;
py::class_<ComponentType> define_pybinding_ComponentType(py::module& m);
} // namespace API
} // namespace MCell

#endif // API_GEN_COMPONENT_TYPE_H
