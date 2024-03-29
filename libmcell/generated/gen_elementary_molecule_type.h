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

#ifndef API_GEN_ELEMENTARY_MOLECULE_TYPE_H
#define API_GEN_ELEMENTARY_MOLECULE_TYPE_H

#include "api/api_common.h"
#include "api/base_data_class.h"

namespace MCell {
namespace API {

class ElementaryMoleculeType;
class Component;
class ComponentType;
class ElementaryMolecule;
class PythonExportContext;

#define ELEMENTARY_MOLECULE_TYPE_CTOR() \
    ElementaryMoleculeType( \
        const std::string& name_, \
        const std::vector<std::shared_ptr<ComponentType>> components_ = std::vector<std::shared_ptr<ComponentType>>(), \
        const double diffusion_constant_2d_ = FLT_UNSET, \
        const double diffusion_constant_3d_ = FLT_UNSET, \
        const double custom_time_step_ = FLT_UNSET, \
        const double custom_space_step_ = FLT_UNSET, \
        const bool target_only_ = false \
    ) { \
      class_name = "ElementaryMoleculeType"; \
      name = name_; \
      components = components_; \
      diffusion_constant_2d = diffusion_constant_2d_; \
      diffusion_constant_3d = diffusion_constant_3d_; \
      custom_time_step = custom_time_step_; \
      custom_space_step = custom_space_step_; \
      target_only = target_only_; \
      postprocess_in_ctor(); \
      check_semantics(); \
    } \
    ElementaryMoleculeType(DefaultCtorArgType) : \
      GenElementaryMoleculeType(DefaultCtorArgType()) { \
      set_all_attributes_as_default_or_unset(); \
      set_all_custom_attributes_to_default(); \
    }

class GenElementaryMoleculeType: public BaseDataClass {
public:
  GenElementaryMoleculeType() {
  }
  GenElementaryMoleculeType(DefaultCtorArgType) {
  }
  void postprocess_in_ctor() override {}
  void check_semantics() const override;
  void set_initialized() override;
  void set_all_attributes_as_default_or_unset() override;

  std::shared_ptr<ElementaryMoleculeType> copy_elementary_molecule_type() const;
  std::shared_ptr<ElementaryMoleculeType> deepcopy_elementary_molecule_type(py::dict = py::dict()) const;
  virtual bool __eq__(const ElementaryMoleculeType& other) const;
  virtual bool eq_nonarray_attributes(const ElementaryMoleculeType& other, const bool ignore_name = false) const;
  bool operator == (const ElementaryMoleculeType& other) const { return __eq__(other);}
  bool operator != (const ElementaryMoleculeType& other) const { return !__eq__(other);}
  std::string to_str(const bool all_details=false, const std::string ind="") const override;

  std::string export_to_python(std::ostream& out, PythonExportContext& ctx) override;
  virtual std::string export_vec_components(std::ostream& out, PythonExportContext& ctx, const std::string& parent_name);


  // --- attributes ---
  std::vector<std::shared_ptr<ComponentType>> components;
  virtual void set_components(const std::vector<std::shared_ptr<ComponentType>> new_components_) {
    if (initialized) {
      throw RuntimeError("Value 'components' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    components = new_components_;
  }
  virtual std::vector<std::shared_ptr<ComponentType>>& get_components() {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return components;
  }

  double diffusion_constant_2d;
  virtual void set_diffusion_constant_2d(const double new_diffusion_constant_2d_) {
    if (initialized) {
      throw RuntimeError("Value 'diffusion_constant_2d' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    diffusion_constant_2d = new_diffusion_constant_2d_;
  }
  virtual double get_diffusion_constant_2d() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return diffusion_constant_2d;
  }

  double diffusion_constant_3d;
  virtual void set_diffusion_constant_3d(const double new_diffusion_constant_3d_) {
    if (initialized) {
      throw RuntimeError("Value 'diffusion_constant_3d' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    diffusion_constant_3d = new_diffusion_constant_3d_;
  }
  virtual double get_diffusion_constant_3d() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return diffusion_constant_3d;
  }

  double custom_time_step;
  virtual void set_custom_time_step(const double new_custom_time_step_) {
    if (initialized) {
      throw RuntimeError("Value 'custom_time_step' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    custom_time_step = new_custom_time_step_;
  }
  virtual double get_custom_time_step() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return custom_time_step;
  }

  double custom_space_step;
  virtual void set_custom_space_step(const double new_custom_space_step_) {
    if (initialized) {
      throw RuntimeError("Value 'custom_space_step' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    custom_space_step = new_custom_space_step_;
  }
  virtual double get_custom_space_step() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return custom_space_step;
  }

  bool target_only;
  virtual void set_target_only(const bool new_target_only_) {
    if (initialized) {
      throw RuntimeError("Value 'target_only' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    target_only = new_target_only_;
  }
  virtual bool get_target_only() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return target_only;
  }

  // --- methods ---
  virtual std::shared_ptr<ElementaryMolecule> inst(const std::vector<std::shared_ptr<Component>> components = std::vector<std::shared_ptr<Component>>(), const std::string& compartment_name = STR_UNSET) = 0;
  virtual std::string to_bngl_str() const = 0;
}; // GenElementaryMoleculeType

class ElementaryMoleculeType;
py::class_<ElementaryMoleculeType> define_pybinding_ElementaryMoleculeType(py::module& m);
} // namespace API
} // namespace MCell

#endif // API_GEN_ELEMENTARY_MOLECULE_TYPE_H
