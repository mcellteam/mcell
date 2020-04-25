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
        const bool accurate_3d_reactions_ = true, \
        const bool center_molecules_on_grid_ = false, \
        const bool micrposcopic_reversibility_ = false \
    ) { \
      class_name = "Config"; \
      time_step = time_step_; \
      surface_grid_density = surface_grid_density_; \
      accurate_3d_reactions = accurate_3d_reactions_; \
      center_molecules_on_grid = center_molecules_on_grid_; \
      micrposcopic_reversibility = micrposcopic_reversibility_; \
      ctor_postprocess();\
    }

class GenConfig: public BaseDataClass {
public:
  void ctor_postprocess() override {}
  SemRes check_semantics(std::ostream& out) const override;
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

  bool accurate_3d_reactions;
  virtual void set_accurate_3d_reactions(const bool new_accurate_3d_reactions_) {
    accurate_3d_reactions = new_accurate_3d_reactions_;
  }
  virtual bool get_accurate_3d_reactions() const {
    return accurate_3d_reactions;
  }

  bool center_molecules_on_grid;
  virtual void set_center_molecules_on_grid(const bool new_center_molecules_on_grid_) {
    center_molecules_on_grid = new_center_molecules_on_grid_;
  }
  virtual bool get_center_molecules_on_grid() const {
    return center_molecules_on_grid;
  }

  bool micrposcopic_reversibility;
  virtual void set_micrposcopic_reversibility(const bool new_micrposcopic_reversibility_) {
    micrposcopic_reversibility = new_micrposcopic_reversibility_;
  }
  virtual bool get_micrposcopic_reversibility() const {
    return micrposcopic_reversibility;
  }

  // --- methods ---
}; // GenConfig

class Config;
py::class_<Config> define_pybinding_Config(py::module& m);
} // namespace API
} // namespace MCell

#endif // API_GEN_CONFIG_H
