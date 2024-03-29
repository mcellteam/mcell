/******************************************************************************
 *
 * Copyright (C) 2021 by
 * The Salk Institute for Biological Studies
 *
 * Use of this source code is governed by an MIT-style
 * license that can be found in the LICENSE file or at
 * https://opensource.org/licenses/MIT.
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
      postprocess_in_ctor(); \
      check_semantics(); \
    } \
    ComponentType(DefaultCtorArgType) : \
      GenComponentType(DefaultCtorArgType()) { \
      set_all_attributes_as_default_or_unset(); \
      set_all_custom_attributes_to_default(); \
    }

class GenComponentType: public BaseDataClass {
public:
  GenComponentType() {
  }
  GenComponentType(DefaultCtorArgType) {
  }
  void postprocess_in_ctor() override {}
  void check_semantics() const override;
  void set_initialized() override;
  void set_all_attributes_as_default_or_unset() override;

  std::shared_ptr<ComponentType> copy_component_type() const;
  std::shared_ptr<ComponentType> deepcopy_component_type(py::dict = py::dict()) const;
  virtual bool __eq__(const ComponentType& other) const;
  virtual bool eq_nonarray_attributes(const ComponentType& other, const bool ignore_name = false) const;
  bool operator == (const ComponentType& other) const { return __eq__(other);}
  bool operator != (const ComponentType& other) const { return !__eq__(other);}
  std::string to_str(const bool all_details=false, const std::string ind="") const override;

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
