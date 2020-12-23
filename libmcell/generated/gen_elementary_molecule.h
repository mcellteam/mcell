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

#ifndef API_GEN_ELEMENTARY_MOLECULE_H
#define API_GEN_ELEMENTARY_MOLECULE_H

#include "api/common.h"
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
        const std::vector<std::shared_ptr<Component>> components_ = std::vector<std::shared_ptr<Component>>() \
    ) { \
      class_name = "ElementaryMolecule"; \
      elementary_molecule_type = elementary_molecule_type_; \
      components = components_; \
      postprocess_in_ctor();\
      check_semantics();\
    }

class GenElementaryMolecule: public BaseDataClass {
public:
  void postprocess_in_ctor() override {}
  void check_semantics() const override;
  void set_initialized() override;
  void set_all_attributes_as_default_or_unset() override;

  virtual bool __eq__(const ElementaryMolecule& other) const;
  virtual bool eq_nonarray_attributes(const ElementaryMolecule& other, const bool ignore_name = false) const;
  bool operator == (const ElementaryMolecule& other) const { return __eq__(other);}
  bool operator != (const ElementaryMolecule& other) const { return !__eq__(other);}
  std::string to_str(const std::string ind="") const override;

  std::string export_to_python(std::ostream& out, PythonExportContext& ctx) const override;
  virtual std::string export_vec_components(std::ostream& out, PythonExportContext& ctx, const std::string& parent_name) const;


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
  virtual std::vector<std::shared_ptr<Component>> get_components() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return components;
  }

  // --- methods ---
  virtual std::string to_bngl_str() const = 0;
}; // GenElementaryMolecule

class ElementaryMolecule;
py::class_<ElementaryMolecule> define_pybinding_ElementaryMolecule(py::module& m);
} // namespace API
} // namespace MCell

#endif // API_GEN_ELEMENTARY_MOLECULE_H
