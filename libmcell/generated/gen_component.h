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

#ifndef API_GEN_COMPONENT_H
#define API_GEN_COMPONENT_H

#include "api/common.h"
#include "api/base_data_class.h"

namespace MCell {
namespace API {

class Component;
class ComponentType;

#define COMPONENT_CTOR() \
    Component( \
        std::shared_ptr<ComponentType> component_type_, \
        const std::string& state_ = "STATE_UNSET", \
        const int bond_ = BOND_UNBOUND \
    ) { \
      class_name = "Component"; \
      component_type = component_type_; \
      state = state_; \
      bond = bond_; \
      postprocess_in_ctor();\
      check_semantics();\
    }

class GenComponent: public BaseDataClass {
public:
  void postprocess_in_ctor() override {}
  void check_semantics() const override;
  void set_initialized() override;
  void set_all_attributes_as_default_or_unset() override;

  virtual bool __eq__(const Component& other) const;
  virtual bool eq_nonarray_attributes(const Component& other, const bool ignore_name = false) const;
  bool operator == (const Component& other) const { return __eq__(other);}
  bool operator != (const Component& other) const { return !__eq__(other);}
  std::string to_str(const std::string ind="") const override;

  std::string export_to_python(std::ostream& out) const override;


  // --- attributes ---
  std::shared_ptr<ComponentType> component_type;
  virtual void set_component_type(std::shared_ptr<ComponentType> new_component_type_) {
    if (initialized) {
      throw RuntimeError("Value 'component_type' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    component_type = new_component_type_;
  }
  virtual std::shared_ptr<ComponentType> get_component_type() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return component_type;
  }

  std::string state;
  virtual void set_state(const std::string& new_state_) {
    if (initialized) {
      throw RuntimeError("Value 'state' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    state = new_state_;
  }
  virtual const std::string& get_state() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return state;
  }

  int bond;
  virtual void set_bond(const int new_bond_) {
    if (initialized) {
      throw RuntimeError("Value 'bond' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    bond = new_bond_;
  }
  virtual int get_bond() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return bond;
  }

  // --- methods ---
  virtual std::string to_bngl_str() const = 0;
}; // GenComponent

class Component;
py::class_<Component> define_pybinding_Component(py::module& m);
} // namespace API
} // namespace MCell

#endif // API_GEN_COMPONENT_H
