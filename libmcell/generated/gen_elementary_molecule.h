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

#ifndef API_GEN_ELEMENTARY_MOLECULE_H
#define API_GEN_ELEMENTARY_MOLECULE_H

#include "api/api_common.h"
#include "api/base_data_class.h"

namespace MCell {
namespace API {

class ElementaryMolecule;
class Component;
class ElementaryMoleculeType;
class PythonExportContext;

#define ELEMENTARY_MOLECULE_CTOR() \
    ElementaryMolecule( \
        std::shared_ptr<ElementaryMoleculeType> elementary_molecule_type_, \
        const std::vector<std::shared_ptr<Component>> components_ = std::vector<std::shared_ptr<Component>>(), \
        const std::string& compartment_name_ = STR_UNSET \
    ) { \
      class_name = "ElementaryMolecule"; \
      elementary_molecule_type = elementary_molecule_type_; \
      components = components_; \
      compartment_name = compartment_name_; \
      postprocess_in_ctor(); \
      check_semantics(); \
    } \
    ElementaryMolecule(DefaultCtorArgType) : \
      GenElementaryMolecule(DefaultCtorArgType()) { \
      set_all_attributes_as_default_or_unset(); \
      set_all_custom_attributes_to_default(); \
    }

class GenElementaryMolecule: public BaseDataClass {
public:
  GenElementaryMolecule() {
  }
  GenElementaryMolecule(DefaultCtorArgType) {
  }
  void postprocess_in_ctor() override {}
  void check_semantics() const override;
  void set_initialized() override;
  void set_all_attributes_as_default_or_unset() override;

  std::shared_ptr<ElementaryMolecule> copy_elementary_molecule() const;
  std::shared_ptr<ElementaryMolecule> deepcopy_elementary_molecule(py::dict = py::dict()) const;
  virtual bool __eq__(const ElementaryMolecule& other) const;
  virtual bool eq_nonarray_attributes(const ElementaryMolecule& other, const bool ignore_name = false) const;
  bool operator == (const ElementaryMolecule& other) const { return __eq__(other);}
  bool operator != (const ElementaryMolecule& other) const { return !__eq__(other);}
  std::string to_str(const bool all_details=false, const std::string ind="") const override;

  std::string export_to_python(std::ostream& out, PythonExportContext& ctx) override;
  virtual std::string export_vec_components(std::ostream& out, PythonExportContext& ctx, const std::string& parent_name);


  // --- attributes ---
  std::shared_ptr<ElementaryMoleculeType> elementary_molecule_type;
  virtual void set_elementary_molecule_type(std::shared_ptr<ElementaryMoleculeType> new_elementary_molecule_type_) {
    if (initialized) {
      throw RuntimeError("Value 'elementary_molecule_type' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    elementary_molecule_type = new_elementary_molecule_type_;
  }
  virtual std::shared_ptr<ElementaryMoleculeType> get_elementary_molecule_type() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return elementary_molecule_type;
  }

  std::vector<std::shared_ptr<Component>> components;
  virtual void set_components(const std::vector<std::shared_ptr<Component>> new_components_) {
    if (initialized) {
      throw RuntimeError("Value 'components' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    components = new_components_;
  }
  virtual std::vector<std::shared_ptr<Component>>& get_components() {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return components;
  }

  std::string compartment_name;
  virtual void set_compartment_name(const std::string& new_compartment_name_) {
    if (initialized) {
      throw RuntimeError("Value 'compartment_name' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    compartment_name = new_compartment_name_;
  }
  virtual const std::string& get_compartment_name() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return compartment_name;
  }

  // --- methods ---
  virtual std::string to_bngl_str(const bool with_compartment = true) const = 0;
}; // GenElementaryMolecule

class ElementaryMolecule;
py::class_<ElementaryMolecule> define_pybinding_ElementaryMolecule(py::module& m);
} // namespace API
} // namespace MCell

#endif // API_GEN_ELEMENTARY_MOLECULE_H
