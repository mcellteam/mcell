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
        const float_t total_iterations_ = 1000000, \
        const bool check_overlapped_walls_ = true, \
        const int reaction_class_cleanup_periodicity_ = 500, \
        const int species_cleanup_periodicity_ = 10000, \
        const bool sort_molecules_ = false, \
        const int memory_limit_gb_ = -1, \
        const uint64_t initial_iteration_ = 0, \
        const float_t initial_time_ = 0, \
        std::shared_ptr<RngState> initial_rng_state_ = nullptr, \
        const bool append_to_count_output_data_ = false \
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
      total_iterations = total_iterations_; \
      check_overlapped_walls = check_overlapped_walls_; \
      reaction_class_cleanup_periodicity = reaction_class_cleanup_periodicity_; \
      species_cleanup_periodicity = species_cleanup_periodicity_; \
      sort_molecules = sort_molecules_; \
      memory_limit_gb = memory_limit_gb_; \
      initial_iteration = initial_iteration_; \
      initial_time = initial_time_; \
      initial_rng_state = initial_rng_state_; \
      append_to_count_output_data = append_to_count_output_data_; \
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

  float_t total_iterations;
  virtual void set_total_iterations(const float_t new_total_iterations_) {
    if (initialized) {
      throw RuntimeError("Value 'total_iterations' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    total_iterations = new_total_iterations_;
  }
  virtual float_t get_total_iterations() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return total_iterations;
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

  int reaction_class_cleanup_periodicity;
  virtual void set_reaction_class_cleanup_periodicity(const int new_reaction_class_cleanup_periodicity_) {
    if (initialized) {
      throw RuntimeError("Value 'reaction_class_cleanup_periodicity' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    reaction_class_cleanup_periodicity = new_reaction_class_cleanup_periodicity_;
  }
  virtual int get_reaction_class_cleanup_periodicity() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return reaction_class_cleanup_periodicity;
  }

  int species_cleanup_periodicity;
  virtual void set_species_cleanup_periodicity(const int new_species_cleanup_periodicity_) {
    if (initialized) {
      throw RuntimeError("Value 'species_cleanup_periodicity' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    species_cleanup_periodicity = new_species_cleanup_periodicity_;
  }
  virtual int get_species_cleanup_periodicity() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return species_cleanup_periodicity;
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

  uint64_t initial_iteration;
  virtual void set_initial_iteration(const uint64_t new_initial_iteration_) {
    if (initialized) {
      throw RuntimeError("Value 'initial_iteration' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    initial_iteration = new_initial_iteration_;
  }
  virtual uint64_t get_initial_iteration() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return initial_iteration;
  }

  float_t initial_time;
  virtual void set_initial_time(const float_t new_initial_time_) {
    if (initialized) {
      throw RuntimeError("Value 'initial_time' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    initial_time = new_initial_time_;
  }
  virtual float_t get_initial_time() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return initial_time;
  }

  std::shared_ptr<RngState> initial_rng_state;
  virtual void set_initial_rng_state(std::shared_ptr<RngState> new_initial_rng_state_) {
    if (initialized) {
      throw RuntimeError("Value 'initial_rng_state' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    initial_rng_state = new_initial_rng_state_;
  }
  virtual std::shared_ptr<RngState> get_initial_rng_state() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return initial_rng_state;
  }

  bool append_to_count_output_data;
  virtual void set_append_to_count_output_data(const bool new_append_to_count_output_data_) {
    if (initialized) {
      throw RuntimeError("Value 'append_to_count_output_data' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    append_to_count_output_data = new_append_to_count_output_data_;
  }
  virtual bool get_append_to_count_output_data() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return append_to_count_output_data;
  }

  // --- methods ---
}; // GenConfig

class Config;
py::class_<Config> define_pybinding_Config(py::module& m);
} // namespace API
} // namespace MCell

#endif // API_GEN_CONFIG_H
