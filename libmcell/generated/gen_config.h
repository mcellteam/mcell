/******************************************************************************
 *
 * Copyright (C) 2020 by
 * The Salk Institute for Biological Studies and
 * Pittsburgh Supercomputing Center, Carnegie Mellon University
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

#ifndef API_GEN_CONFIG_H
#define API_GEN_CONFIG_H

#include "../api/common.h"

namespace MCell {
namespace API {

#define CONFIG_CTOR() \
    Config( \
        const float_t time_step_ = 1e-6, \
        const float_t surface_grid_density_ = 10000, \
        const bool center_molecules_on_grid_ = false, \
        const bool microscopic_reversibility_ = false \
    ) { \
      class_name = "Config"; \
      time_step = time_step_; \
      surface_grid_density = surface_grid_density_; \
      center_molecules_on_grid = center_molecules_on_grid_; \
      microscopic_reversibility = microscopic_reversibility_; \
      postprocess_in_ctor();\
      check_semantics();\
    }

class GenConfig: public BaseDataClass {
public:
  void postprocess_in_ctor() override {}
  void check_semantics() const override;
  std::string to_str(const std::string ind="") const override;

  // --- attributes ---
  float_t time_step;
  virtual void set_time_step(const float_t new_time_step_) {
    time_step = new_time_step_;
  }
  virtual float_t get_time_step() const {
    return time_step;
  }

  float_t surface_grid_density;
  virtual void set_surface_grid_density(const float_t new_surface_grid_density_) {
    surface_grid_density = new_surface_grid_density_;
  }
  virtual float_t get_surface_grid_density() const {
    return surface_grid_density;
  }

  bool center_molecules_on_grid;
  virtual void set_center_molecules_on_grid(const bool new_center_molecules_on_grid_) {
    center_molecules_on_grid = new_center_molecules_on_grid_;
  }
  virtual bool get_center_molecules_on_grid() const {
    return center_molecules_on_grid;
  }

  bool microscopic_reversibility;
  virtual void set_microscopic_reversibility(const bool new_microscopic_reversibility_) {
    microscopic_reversibility = new_microscopic_reversibility_;
  }
  virtual bool get_microscopic_reversibility() const {
    return microscopic_reversibility;
  }

  // --- methods ---
}; // GenConfig

class Config;
py::class_<Config> define_pybinding_Config(py::module& m);
} // namespace API
} // namespace MCell

#endif // API_GEN_CONFIG_H
