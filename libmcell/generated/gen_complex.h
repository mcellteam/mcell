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

#ifndef API_GEN_COMPLEX_H
#define API_GEN_COMPLEX_H

#include "api/common.h"
#include "api/base_data_class.h"

namespace MCell {
namespace API {

class Complex;
class ElementaryMolecule;
class Species;
class PythonExportContext;

#define COMPLEX_CTOR() \
    Complex( \
        const std::string& name_ = STR_UNSET, \
        const std::vector<std::shared_ptr<ElementaryMolecule>> elementary_molecules_ = std::vector<std::shared_ptr<ElementaryMolecule>>(), \
        const Orientation orientation_ = Orientation::DEFAULT, \
        const std::string& compartment_name_ = STR_UNSET \
    ) { \
      class_name = "Complex"; \
      name = name_; \
      elementary_molecules = elementary_molecules_; \
      orientation = orientation_; \
      compartment_name = compartment_name_; \
      postprocess_in_ctor();\
      check_semantics();\
    }

class GenComplex: public BaseDataClass {
public:
  void postprocess_in_ctor() override {}
  void check_semantics() const override;
  void set_initialized() override;
  void set_all_attributes_as_default_or_unset() override;

  virtual bool __eq__(const Complex& other) const;
  virtual bool eq_nonarray_attributes(const Complex& other, const bool ignore_name = false) const;
  bool operator == (const Complex& other) const { return __eq__(other);}
  bool operator != (const Complex& other) const { return !__eq__(other);}
  std::string to_str(const std::string ind="") const override;

  std::string export_to_python(std::ostream& out, PythonExportContext& ctx) override;
  virtual std::string export_vec_elementary_molecules(std::ostream& out, PythonExportContext& ctx, const std::string& parent_name);


  // --- attributes ---
  std::vector<std::shared_ptr<ElementaryMolecule>> elementary_molecules;
  virtual void set_elementary_molecules(const std::vector<std::shared_ptr<ElementaryMolecule>> new_elementary_molecules_) {
    if (initialized) {
      throw RuntimeError("Value 'elementary_molecules' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    elementary_molecules = new_elementary_molecules_;
  }
  virtual std::vector<std::shared_ptr<ElementaryMolecule>> get_elementary_molecules() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return elementary_molecules;
  }

  Orientation orientation;
  virtual void set_orientation(const Orientation new_orientation_) {
    if (initialized) {
      throw RuntimeError("Value 'orientation' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    orientation = new_orientation_;
  }
  virtual Orientation get_orientation() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return orientation;
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
  virtual std::string to_bngl_str() const = 0;
  virtual std::shared_ptr<Species> as_species() = 0;
}; // GenComplex

class Complex;
py::class_<Complex> define_pybinding_Complex(py::module& m);
} // namespace API
} // namespace MCell

#endif // API_GEN_COMPLEX_H
