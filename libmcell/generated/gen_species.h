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

#ifndef API_GEN_SPECIES_H
#define API_GEN_SPECIES_H

#include "api/api_common.h"
#include "api/complex.h"


namespace MCell {
namespace API {

class Complex;
class ElementaryMolecule;
class Species;
class PythonExportContext;

#define SPECIES_CTOR() \
    Species( \
        const std::string& name_ = STR_UNSET, \
        const double diffusion_constant_2d_ = FLT_UNSET, \
        const double diffusion_constant_3d_ = FLT_UNSET, \
        const double custom_time_step_ = FLT_UNSET, \
        const double custom_space_step_ = FLT_UNSET, \
        const bool target_only_ = false, \
        const std::vector<std::shared_ptr<ElementaryMolecule>> elementary_molecules_ = std::vector<std::shared_ptr<ElementaryMolecule>>(), \
        const Orientation orientation_ = Orientation::DEFAULT, \
        const std::string& compartment_name_ = STR_UNSET \
    )  : GenSpecies(name_,elementary_molecules_,orientation_,compartment_name_) { \
      class_name = "Species"; \
      name = name_; \
      diffusion_constant_2d = diffusion_constant_2d_; \
      diffusion_constant_3d = diffusion_constant_3d_; \
      custom_time_step = custom_time_step_; \
      custom_space_step = custom_space_step_; \
      target_only = target_only_; \
      elementary_molecules = elementary_molecules_; \
      orientation = orientation_; \
      compartment_name = compartment_name_; \
      postprocess_in_ctor(); \
      check_semantics(); \
    } \
    Species(DefaultCtorArgType) : \
      GenSpecies(DefaultCtorArgType()) { \
      set_all_attributes_as_default_or_unset(); \
      set_all_custom_attributes_to_default(); \
    }

class GenSpecies: public Complex {
public:
  GenSpecies( 
      const std::string& name_ = STR_UNSET, 
      const std::vector<std::shared_ptr<ElementaryMolecule>> elementary_molecules_ = std::vector<std::shared_ptr<ElementaryMolecule>>(), 
      const Orientation orientation_ = Orientation::DEFAULT, 
      const std::string& compartment_name_ = STR_UNSET 
  )  : Complex(name_,elementary_molecules_,orientation_,compartment_name_)  {
  }
  GenSpecies(DefaultCtorArgType) {
  }
  void postprocess_in_ctor() override {}
  void check_semantics() const override;
  void set_initialized() override;
  void set_all_attributes_as_default_or_unset() override;

  std::shared_ptr<Species> copy_species() const;
  std::shared_ptr<Species> deepcopy_species(py::dict = py::dict()) const;
  virtual bool __eq__(const Species& other) const;
  virtual bool eq_nonarray_attributes(const Species& other, const bool ignore_name = false) const;
  bool operator == (const Species& other) const { return __eq__(other);}
  bool operator != (const Species& other) const { return !__eq__(other);}
  std::string to_str(const bool all_details=false, const std::string ind="") const override;

  virtual std::string export_to_python(std::ostream& out, PythonExportContext& ctx);
  virtual std::string export_vec_elementary_molecules(std::ostream& out, PythonExportContext& ctx, const std::string& parent_name);


  // --- attributes ---
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
  virtual std::shared_ptr<Complex> inst(const Orientation orientation = Orientation::DEFAULT, const std::string& compartment_name = STR_UNSET) = 0;
}; // GenSpecies

class Species;
py::class_<Species> define_pybinding_Species(py::module& m);
} // namespace API
} // namespace MCell

#endif // API_GEN_SPECIES_H
