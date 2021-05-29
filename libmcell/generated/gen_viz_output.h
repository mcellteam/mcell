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

#ifndef API_GEN_VIZ_OUTPUT_H
#define API_GEN_VIZ_OUTPUT_H

#include "api/api_common.h"
#include "api/base_data_class.h"

namespace MCell {
namespace API {

class VizOutput;
class Species;
class PythonExportContext;

#define VIZ_OUTPUT_CTOR() \
    VizOutput( \
        const std::string& output_files_prefix_, \
        const std::vector<std::shared_ptr<Species>> species_list_ = std::vector<std::shared_ptr<Species>>(), \
        const VizMode mode_ = VizMode::ASCII, \
        const double every_n_timesteps_ = 1 \
    ) { \
      class_name = "VizOutput"; \
      output_files_prefix = output_files_prefix_; \
      species_list = species_list_; \
      mode = mode_; \
      every_n_timesteps = every_n_timesteps_; \
      postprocess_in_ctor(); \
      check_semantics(); \
    } \
    VizOutput(DefaultCtorArgType) : \
      GenVizOutput(DefaultCtorArgType()) { \
      set_all_attributes_as_default_or_unset(); \
    }

class GenVizOutput: public BaseDataClass {
public:
  GenVizOutput() {
  }
  GenVizOutput(DefaultCtorArgType) {
  }
  void postprocess_in_ctor() override {}
  void check_semantics() const override;
  void set_initialized() override;
  void set_all_attributes_as_default_or_unset() override;

  std::shared_ptr<VizOutput> copy_viz_output() const;
  std::shared_ptr<VizOutput> deepcopy_viz_output(py::dict = py::dict()) const;
  virtual bool __eq__(const VizOutput& other) const;
  virtual bool eq_nonarray_attributes(const VizOutput& other, const bool ignore_name = false) const;
  bool operator == (const VizOutput& other) const { return __eq__(other);}
  bool operator != (const VizOutput& other) const { return !__eq__(other);}
  std::string to_str(const bool all_details=false, const std::string ind="") const override;

  std::string export_to_python(std::ostream& out, PythonExportContext& ctx) override;
  virtual std::string export_vec_species_list(std::ostream& out, PythonExportContext& ctx, const std::string& parent_name);


  // --- attributes ---
  std::string output_files_prefix;
  virtual void set_output_files_prefix(const std::string& new_output_files_prefix_) {
    if (initialized) {
      throw RuntimeError("Value 'output_files_prefix' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    output_files_prefix = new_output_files_prefix_;
  }
  virtual const std::string& get_output_files_prefix() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return output_files_prefix;
  }

  std::vector<std::shared_ptr<Species>> species_list;
  virtual void set_species_list(const std::vector<std::shared_ptr<Species>> new_species_list_) {
    if (initialized) {
      throw RuntimeError("Value 'species_list' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    species_list = new_species_list_;
  }
  virtual std::vector<std::shared_ptr<Species>>& get_species_list() {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return species_list;
  }

  VizMode mode;
  virtual void set_mode(const VizMode new_mode_) {
    if (initialized) {
      throw RuntimeError("Value 'mode' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    mode = new_mode_;
  }
  virtual VizMode get_mode() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return mode;
  }

  double every_n_timesteps;
  virtual void set_every_n_timesteps(const double new_every_n_timesteps_) {
    if (initialized) {
      throw RuntimeError("Value 'every_n_timesteps' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    every_n_timesteps = new_every_n_timesteps_;
  }
  virtual double get_every_n_timesteps() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return every_n_timesteps;
  }

  // --- methods ---
}; // GenVizOutput

class VizOutput;
py::class_<VizOutput> define_pybinding_VizOutput(py::module& m);
} // namespace API
} // namespace MCell

#endif // API_GEN_VIZ_OUTPUT_H
