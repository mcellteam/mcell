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

#ifndef API_GEN_WARNINGS_H
#define API_GEN_WARNINGS_H

#include "api/api_common.h"
#include "api/base_data_class.h"

namespace MCell {
namespace API {

class Warnings;
class PythonExportContext;

#define WARNINGS_CTOR() \
    Warnings( \
        const WarningLevel high_reaction_probability_ = WarningLevel::IGNORE, \
        const WarningLevel molecule_placement_failure_ = WarningLevel::ERROR \
    ) { \
      class_name = "Warnings"; \
      high_reaction_probability = high_reaction_probability_; \
      molecule_placement_failure = molecule_placement_failure_; \
      postprocess_in_ctor(); \
      check_semantics(); \
    } \
    Warnings(DefaultCtorArgType) : \
      GenWarnings(DefaultCtorArgType()) { \
      set_all_attributes_as_default_or_unset(); \
      set_all_custom_attributes_to_default(); \
    }

class GenWarnings: public BaseDataClass {
public:
  GenWarnings() {
  }
  GenWarnings(DefaultCtorArgType) {
  }
  void postprocess_in_ctor() override {}
  void check_semantics() const override;
  void set_initialized() override;
  void set_all_attributes_as_default_or_unset() override;

  std::shared_ptr<Warnings> copy_warnings() const;
  std::shared_ptr<Warnings> deepcopy_warnings(py::dict = py::dict()) const;
  virtual bool __eq__(const Warnings& other) const;
  virtual bool eq_nonarray_attributes(const Warnings& other, const bool ignore_name = false) const;
  bool operator == (const Warnings& other) const { return __eq__(other);}
  bool operator != (const Warnings& other) const { return !__eq__(other);}
  std::string to_str(const bool all_details=false, const std::string ind="") const override;

  std::string export_to_python(std::ostream& out, PythonExportContext& ctx) override;


  // --- attributes ---
  WarningLevel high_reaction_probability;
  virtual void set_high_reaction_probability(const WarningLevel new_high_reaction_probability_) {
    if (initialized) {
      throw RuntimeError("Value 'high_reaction_probability' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    high_reaction_probability = new_high_reaction_probability_;
  }
  virtual WarningLevel get_high_reaction_probability() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return high_reaction_probability;
  }

  WarningLevel molecule_placement_failure;
  virtual void set_molecule_placement_failure(const WarningLevel new_molecule_placement_failure_) {
    if (initialized) {
      throw RuntimeError("Value 'molecule_placement_failure' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    molecule_placement_failure = new_molecule_placement_failure_;
  }
  virtual WarningLevel get_molecule_placement_failure() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return molecule_placement_failure;
  }

  // --- methods ---
}; // GenWarnings

class Warnings;
py::class_<Warnings> define_pybinding_Warnings(py::module& m);
} // namespace API
} // namespace MCell

#endif // API_GEN_WARNINGS_H
