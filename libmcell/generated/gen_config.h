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

#ifndef API_GEN_CONFIG_H
#define API_GEN_CONFIG_H

#include "../api/common.h"
#include "../api/base_data_class.h"

namespace MCell {
namespace API {

#define CONFIG_CTOR() \
    Config( \
        const int seed_ = 1, \
        const float_t time_step_ = 1e-6, \
        const float_t surface_grid_density_ = 10000, \
        const bool center_molecules_on_grid_ = false, \
        const bool microscopic_reversibility_ = false, \
        const float_t partition_dimension_ = 10, \
        const float_t subpartition_dimension_ = 0.5 \
    ) { \
      class_name = "Config"; \
      seed = seed_; \
      time_step = time_step_; \
      surface_grid_density = surface_grid_density_; \
      center_molecules_on_grid = center_molecules_on_grid_; \
      microscopic_reversibility = microscopic_reversibility_; \
      partition_dimension = partition_dimension_; \
      subpartition_dimension = subpartition_dimension_; \
      postprocess_in_ctor();\
      check_semantics();\
    }

class GenConfig: public BaseDataClass {
public:
  void postprocess_in_ctor() override {}
  void check_semantics() const override;
  std::string to_str(const std::string ind="") const override;

  bool __eq__(const GenConfig& other) const;
  // --- attributes ---
  int seed;
  virtual void set_seed(const int new_seed_) {
    seed = new_seed_;
  }
  virtual int get_seed() const {
    return seed;
  }

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

  float_t partition_dimension;
  virtual void set_partition_dimension(const float_t new_partition_dimension_) {
    partition_dimension = new_partition_dimension_;
  }
  virtual float_t get_partition_dimension() const {
    return partition_dimension;
  }

  float_t subpartition_dimension;
  virtual void set_subpartition_dimension(const float_t new_subpartition_dimension_) {
    subpartition_dimension = new_subpartition_dimension_;
  }
  virtual float_t get_subpartition_dimension() const {
    return subpartition_dimension;
  }

  // --- methods ---
}; // GenConfig

class Config;
py::class_<Config> define_pybinding_Config(py::module& m);
} // namespace API
} // namespace MCell

#endif // API_GEN_CONFIG_H
