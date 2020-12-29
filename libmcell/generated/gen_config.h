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

#include "api/common.h"
#include "api/base_data_class.h"

namespace MCell {
namespace API {

class Config;
class RngState;
class PythonExportContext;

#define CONFIG_CTOR() \
    Config( \
        const int seed_ = 1, \
        const float_t time_step_ = 1e-6, \
        const float_t surface_grid_density_ = 10000, \
        const float_t interaction_radius_ = FLT_UNSET, \
        const float_t intermembrane_interaction_radius_ = FLT_UNSET, \
        const float_t vacancy_search_distance_ = 10, \
        const bool center_molecules_on_grid_ = false, \
        const std::vector<float_t> initial_partition_origin_ = std::vector<float_t>(), \
        const float_t partition_dimension_ = 10, \
        const float_t subpartition_dimension_ = 0.5, \
        const float_t total_iterations_hint_ = 1000000, \
        const bool check_overlapped_walls_ = true, \
        const bool sort_molecules_ = false, \
        const int memory_limit_gb_ = -1, \
        std::shared_ptr<RngState> rng_state_ = nullptr \
    ) { \
      class_name = "Config"; \
      seed = seed_; \
      time_step = time_step_; \
      surface_grid_density = surface_grid_density_; \
      interaction_radius = interaction_radius_; \
      intermembrane_interaction_radius = intermembrane_interaction_radius_; \
      vacancy_search_distance = vacancy_search_distance_; \
      center_molecules_on_grid = center_molecules_on_grid_; \
      initial_partition_origin = initial_partition_origin_; \
      partition_dimension = partition_dimension_; \
      subpartition_dimension = subpartition_dimension_; \
      total_iterations_hint = total_iterations_hint_; \
      check_overlapped_walls = check_overlapped_walls_; \
      sort_molecules = sort_molecules_; \
      memory_limit_gb = memory_limit_gb_; \
      rng_state = rng_state_; \
      postprocess_in_ctor();\
      check_semantics();\
    }

class GenConfig: public BaseDataClass {
public:
  void postprocess_in_ctor() override {}
  void check_semantics() const override;
  void set_initialized() override;
  void set_all_attributes_as_default_or_unset() override;

  virtual bool __eq__(const Config& other) const;
  virtual bool eq_nonarray_attributes(const Config& other, const bool ignore_name = false) const;
  bool operator == (const Config& other) const { return __eq__(other);}
  bool operator != (const Config& other) const { return !__eq__(other);}
  std::string to_str(const std::string ind="") const override;

  std::string export_to_python(std::ostream& out, PythonExportContext& ctx) const override;
  virtual std::string export_vec_initial_partition_origin(std::ostream& out, PythonExportContext& ctx, const std::string& parent_name) const;


  // --- attributes ---
  int seed;
  virtual void set_seed(const int new_seed_) {
    if (initialized) {
      throw RuntimeError("Value 'seed' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    seed = new_seed_;
  }
  virtual int get_seed() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return seed;
  }

  float_t time_step;
  virtual void set_time_step(const float_t new_time_step_) {
    if (initialized) {
      throw RuntimeError("Value 'time_step' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    time_step = new_time_step_;
  }
  virtual float_t get_time_step() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return time_step;
  }

  float_t surface_grid_density;
  virtual void set_surface_grid_density(const float_t new_surface_grid_density_) {
    if (initialized) {
      throw RuntimeError("Value 'surface_grid_density' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    surface_grid_density = new_surface_grid_density_;
  }
  virtual float_t get_surface_grid_density() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return surface_grid_density;
  }

  float_t interaction_radius;
  virtual void set_interaction_radius(const float_t new_interaction_radius_) {
    if (initialized) {
      throw RuntimeError("Value 'interaction_radius' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    interaction_radius = new_interaction_radius_;
  }
  virtual float_t get_interaction_radius() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return interaction_radius;
  }

  float_t intermembrane_interaction_radius;
  virtual void set_intermembrane_interaction_radius(const float_t new_intermembrane_interaction_radius_) {
    if (initialized) {
      throw RuntimeError("Value 'intermembrane_interaction_radius' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    intermembrane_interaction_radius = new_intermembrane_interaction_radius_;
  }
  virtual float_t get_intermembrane_interaction_radius() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return intermembrane_interaction_radius;
  }

  float_t vacancy_search_distance;
  virtual void set_vacancy_search_distance(const float_t new_vacancy_search_distance_) {
    if (initialized) {
      throw RuntimeError("Value 'vacancy_search_distance' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    vacancy_search_distance = new_vacancy_search_distance_;
  }
  virtual float_t get_vacancy_search_distance() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return vacancy_search_distance;
  }

  bool center_molecules_on_grid;
  virtual void set_center_molecules_on_grid(const bool new_center_molecules_on_grid_) {
    if (initialized) {
      throw RuntimeError("Value 'center_molecules_on_grid' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    center_molecules_on_grid = new_center_molecules_on_grid_;
  }
  virtual bool get_center_molecules_on_grid() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return center_molecules_on_grid;
  }

  std::vector<float_t> initial_partition_origin;
  virtual void set_initial_partition_origin(const std::vector<float_t> new_initial_partition_origin_) {
    if (initialized) {
      throw RuntimeError("Value 'initial_partition_origin' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    initial_partition_origin = new_initial_partition_origin_;
  }
  virtual std::vector<float_t> get_initial_partition_origin() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return initial_partition_origin;
  }

  float_t partition_dimension;
  virtual void set_partition_dimension(const float_t new_partition_dimension_) {
    if (initialized) {
      throw RuntimeError("Value 'partition_dimension' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    partition_dimension = new_partition_dimension_;
  }
  virtual float_t get_partition_dimension() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return partition_dimension;
  }

  float_t subpartition_dimension;
  virtual void set_subpartition_dimension(const float_t new_subpartition_dimension_) {
    if (initialized) {
      throw RuntimeError("Value 'subpartition_dimension' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    subpartition_dimension = new_subpartition_dimension_;
  }
  virtual float_t get_subpartition_dimension() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return subpartition_dimension;
  }

  float_t total_iterations_hint;
  virtual void set_total_iterations_hint(const float_t new_total_iterations_hint_) {
    if (initialized) {
      throw RuntimeError("Value 'total_iterations_hint' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    total_iterations_hint = new_total_iterations_hint_;
  }
  virtual float_t get_total_iterations_hint() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return total_iterations_hint;
  }

  bool check_overlapped_walls;
  virtual void set_check_overlapped_walls(const bool new_check_overlapped_walls_) {
    if (initialized) {
      throw RuntimeError("Value 'check_overlapped_walls' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    check_overlapped_walls = new_check_overlapped_walls_;
  }
  virtual bool get_check_overlapped_walls() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return check_overlapped_walls;
  }

  bool sort_molecules;
  virtual void set_sort_molecules(const bool new_sort_molecules_) {
    if (initialized) {
      throw RuntimeError("Value 'sort_molecules' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    sort_molecules = new_sort_molecules_;
  }
  virtual bool get_sort_molecules() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return sort_molecules;
  }

  int memory_limit_gb;
  virtual void set_memory_limit_gb(const int new_memory_limit_gb_) {
    if (initialized) {
      throw RuntimeError("Value 'memory_limit_gb' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    memory_limit_gb = new_memory_limit_gb_;
  }
  virtual int get_memory_limit_gb() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return memory_limit_gb;
  }

  std::shared_ptr<RngState> rng_state;
  virtual void set_rng_state(std::shared_ptr<RngState> new_rng_state_) {
    if (initialized) {
      throw RuntimeError("Value 'rng_state' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    rng_state = new_rng_state_;
  }
  virtual std::shared_ptr<RngState> get_rng_state() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return rng_state;
  }

  // --- methods ---
}; // GenConfig

class Config;
py::class_<Config> define_pybinding_Config(py::module& m);
} // namespace API
} // namespace MCell

#endif // API_GEN_CONFIG_H
