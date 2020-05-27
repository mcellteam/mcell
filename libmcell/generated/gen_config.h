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
        const float_t interaction_radius_ = FLT_UNSET, \
        const float_t vacancy_search_distance_ = 1e+6, \
        const bool center_molecules_on_grid_ = false, \
        const float_t partition_dimension_ = 10, \
        const float_t subpartition_dimension_ = 0.5, \
        const long total_iterations_hint_ = 1000000 \
    ) { \
      class_name = "Config"; \
      seed = seed_; \
      time_step = time_step_; \
      surface_grid_density = surface_grid_density_; \
      interaction_radius = interaction_radius_; \
      vacancy_search_distance = vacancy_search_distance_; \
      center_molecules_on_grid = center_molecules_on_grid_; \
      partition_dimension = partition_dimension_; \
      subpartition_dimension = subpartition_dimension_; \
      total_iterations_hint = total_iterations_hint_; \
      postprocess_in_ctor();\
      check_semantics();\
    }

class GenConfig: public BaseDataClass {
public:
  void postprocess_in_ctor() override {}
  void check_semantics() const override;
  void set_initialized() override;

  bool __eq__(const GenConfig& other) const;
  std::string to_str(const std::string ind="") const override;

  // --- attributes ---
  int seed;
  virtual void set_seed(const int new_seed_) {
    if (initialized) {
      throw RuntimeError("Value 'seed' of object with name " + name + " (class " + class_name + ")"
                         "cannot be set after model was initialized.");
    }
    seed = new_seed_;
  }
  virtual int get_seed() const {
    return seed;
  }

  float_t time_step;
  virtual void set_time_step(const float_t new_time_step_) {
    if (initialized) {
      throw RuntimeError("Value 'time_step' of object with name " + name + " (class " + class_name + ")"
                         "cannot be set after model was initialized.");
    }
    time_step = new_time_step_;
  }
  virtual float_t get_time_step() const {
    return time_step;
  }

  float_t surface_grid_density;
  virtual void set_surface_grid_density(const float_t new_surface_grid_density_) {
    if (initialized) {
      throw RuntimeError("Value 'surface_grid_density' of object with name " + name + " (class " + class_name + ")"
                         "cannot be set after model was initialized.");
    }
    surface_grid_density = new_surface_grid_density_;
  }
  virtual float_t get_surface_grid_density() const {
    return surface_grid_density;
  }

  float_t interaction_radius;
  virtual void set_interaction_radius(const float_t new_interaction_radius_) {
    if (initialized) {
      throw RuntimeError("Value 'interaction_radius' of object with name " + name + " (class " + class_name + ")"
                         "cannot be set after model was initialized.");
    }
    interaction_radius = new_interaction_radius_;
  }
  virtual float_t get_interaction_radius() const {
    return interaction_radius;
  }

  float_t vacancy_search_distance;
  virtual void set_vacancy_search_distance(const float_t new_vacancy_search_distance_) {
    if (initialized) {
      throw RuntimeError("Value 'vacancy_search_distance' of object with name " + name + " (class " + class_name + ")"
                         "cannot be set after model was initialized.");
    }
    vacancy_search_distance = new_vacancy_search_distance_;
  }
  virtual float_t get_vacancy_search_distance() const {
    return vacancy_search_distance;
  }

  bool center_molecules_on_grid;
  virtual void set_center_molecules_on_grid(const bool new_center_molecules_on_grid_) {
    if (initialized) {
      throw RuntimeError("Value 'center_molecules_on_grid' of object with name " + name + " (class " + class_name + ")"
                         "cannot be set after model was initialized.");
    }
    center_molecules_on_grid = new_center_molecules_on_grid_;
  }
  virtual bool get_center_molecules_on_grid() const {
    return center_molecules_on_grid;
  }

  float_t partition_dimension;
  virtual void set_partition_dimension(const float_t new_partition_dimension_) {
    if (initialized) {
      throw RuntimeError("Value 'partition_dimension' of object with name " + name + " (class " + class_name + ")"
                         "cannot be set after model was initialized.");
    }
    partition_dimension = new_partition_dimension_;
  }
  virtual float_t get_partition_dimension() const {
    return partition_dimension;
  }

  float_t subpartition_dimension;
  virtual void set_subpartition_dimension(const float_t new_subpartition_dimension_) {
    if (initialized) {
      throw RuntimeError("Value 'subpartition_dimension' of object with name " + name + " (class " + class_name + ")"
                         "cannot be set after model was initialized.");
    }
    subpartition_dimension = new_subpartition_dimension_;
  }
  virtual float_t get_subpartition_dimension() const {
    return subpartition_dimension;
  }

  long total_iterations_hint;
  virtual void set_total_iterations_hint(const long new_total_iterations_hint_) {
    if (initialized) {
      throw RuntimeError("Value 'total_iterations_hint' of object with name " + name + " (class " + class_name + ")"
                         "cannot be set after model was initialized.");
    }
    total_iterations_hint = new_total_iterations_hint_;
  }
  virtual long get_total_iterations_hint() const {
    return total_iterations_hint;
  }

  // --- methods ---
}; // GenConfig

class Config;
py::class_<Config> define_pybinding_Config(py::module& m);
} // namespace API
} // namespace MCell

#endif // API_GEN_CONFIG_H
