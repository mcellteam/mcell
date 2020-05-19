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

#include "../api/common.h"
#include "../api/base_data_class.h"

namespace MCell {
namespace API {

class Species;

#define VIZ_OUTPUT_CTOR() \
    VizOutput( \
        const std::string& filename_, \
        const std::vector<std::shared_ptr<Species>> species_list_ = std::vector<std::shared_ptr<Species>>(), \
        const VizMode mode_ = VizMode::Ascii, \
        const int every_n_timesteps_ = 1 \
    ) { \
      class_name = "VizOutput"; \
      filename = filename_; \
      species_list = species_list_; \
      mode = mode_; \
      every_n_timesteps = every_n_timesteps_; \
      postprocess_in_ctor();\
      check_semantics();\
    }

class GenVizOutput: public BaseDataClass {
public:
  void postprocess_in_ctor() override {}
  void check_semantics() const override;
  std::string to_str(const std::string ind="") const override;

  void set_initialized() override;

  bool __eq__(const GenVizOutput& other) const;
  // --- attributes ---
  std::string filename;
  virtual void set_filename(const std::string& new_filename_) {
    if (initialized) {
      throw RuntimeError("Value 'filename' of object with name " + name + " (class " + class_name + ")"
                         "cannot be set after model was initialized.");
    }
    filename = new_filename_;
  }
  virtual const std::string& get_filename() const {
    return filename;
  }

  std::vector<std::shared_ptr<Species>> species_list;
  virtual void set_species_list(const std::vector<std::shared_ptr<Species>> new_species_list_) {
    if (initialized) {
      throw RuntimeError("Value 'species_list' of object with name " + name + " (class " + class_name + ")"
                         "cannot be set after model was initialized.");
    }
    species_list = new_species_list_;
  }
  virtual std::vector<std::shared_ptr<Species>> get_species_list() const {
    return species_list;
  }

  VizMode mode;
  virtual void set_mode(const VizMode new_mode_) {
    if (initialized) {
      throw RuntimeError("Value 'mode' of object with name " + name + " (class " + class_name + ")"
                         "cannot be set after model was initialized.");
    }
    mode = new_mode_;
  }
  virtual VizMode get_mode() const {
    return mode;
  }

  int every_n_timesteps;
  virtual void set_every_n_timesteps(const int new_every_n_timesteps_) {
    if (initialized) {
      throw RuntimeError("Value 'every_n_timesteps' of object with name " + name + " (class " + class_name + ")"
                         "cannot be set after model was initialized.");
    }
    every_n_timesteps = new_every_n_timesteps_;
  }
  virtual int get_every_n_timesteps() const {
    return every_n_timesteps;
  }

  // --- methods ---
}; // GenVizOutput

class VizOutput;
py::class_<VizOutput> define_pybinding_VizOutput(py::module& m);
} // namespace API
} // namespace MCell

#endif // API_GEN_VIZ_OUTPUT_H
