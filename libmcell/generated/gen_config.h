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

#ifndef API_GEN_CONFIG_H
#define API_GEN_CONFIG_H

#include "api/api_common.h"
#include "api/base_data_class.h"

namespace MCell {
namespace API {

class Config;
class RngState;
class PythonExportContext;

#define CONFIG_CTOR() \
    Config( \
        const int seed_ = 1, \
        const double time_step_ = 1e-6, \
        const double surface_grid_density_ = 10000, \
        const double interaction_radius_ = FLT_UNSET, \
        const double intermembrane_interaction_radius_ = FLT_UNSET, \
        const double vacancy_search_distance_ = 10, \
        const bool center_molecules_on_grid_ = false, \
        const double partition_dimension_ = 10, \
        const std::vector<double> initial_partition_origin_ = std::vector<double>(), \
        const double subpartition_dimension_ = 0.5, \
        const double total_iterations_ = 1000000, \
        const bool check_overlapped_walls_ = true, \
        const int reaction_class_cleanup_periodicity_ = 500, \
        const int species_cleanup_periodicity_ = 10000, \
        const int molecules_order_random_shuffle_periodicity_ = 10000, \
        const bool sort_molecules_ = false, \
        const int memory_limit_gb_ = -1, \
        const uint64_t initial_iteration_ = 0, \
        const double initial_time_ = 0, \
        std::shared_ptr<RngState> initial_rng_state_ = nullptr, \
        const bool append_to_count_output_data_ = false, \
        const bool continue_after_sigalrm_ = false \
    ) { \
      class_name = "Config"; \
      seed = seed_; \
      time_step = time_step_; \
      surface_grid_density = surface_grid_density_; \
      interaction_radius = interaction_radius_; \
      intermembrane_interaction_radius = intermembrane_interaction_radius_; \
      vacancy_search_distance = vacancy_search_distance_; \
      center_molecules_on_grid = center_molecules_on_grid_; \
      partition_dimension = partition_dimension_; \
      initial_partition_origin = initial_partition_origin_; \
      subpartition_dimension = subpartition_dimension_; \
      total_iterations = total_iterations_; \
      check_overlapped_walls = check_overlapped_walls_; \
      reaction_class_cleanup_periodicity = reaction_class_cleanup_periodicity_; \
      species_cleanup_periodicity = species_cleanup_periodicity_; \
      molecules_order_random_shuffle_periodicity = molecules_order_random_shuffle_periodicity_; \
      sort_molecules = sort_molecules_; \
      memory_limit_gb = memory_limit_gb_; \
      initial_iteration = initial_iteration_; \
      initial_time = initial_time_; \
      initial_rng_state = initial_rng_state_; \
      append_to_count_output_data = append_to_count_output_data_; \
      continue_after_sigalrm = continue_after_sigalrm_; \
      postprocess_in_ctor(); \
      check_semantics(); \
    } \
    Config(DefaultCtorArgType) : \
      GenConfig(DefaultCtorArgType()) { \
      set_all_attributes_as_default_or_unset(); \
      set_all_custom_attributes_to_default(); \
    }

class GenConfig: public BaseDataClass {
public:
  GenConfig() {
  }
  GenConfig(DefaultCtorArgType) {
  }
  void postprocess_in_ctor() override {}
  void check_semantics() const override;
  void set_initialized() override;
  void set_all_attributes_as_default_or_unset() override;

  std::shared_ptr<Config> copy_config() const;
  std::shared_ptr<Config> deepcopy_config(py::dict = py::dict()) const;
  virtual bool __eq__(const Config& other) const;
  virtual bool eq_nonarray_attributes(const Config& other, const bool ignore_name = false) const;
  bool operator == (const Config& other) const { return __eq__(other);}
  bool operator != (const Config& other) const { return !__eq__(other);}
  std::string to_str(const bool all_details=false, const std::string ind="") const override;

  std::string export_to_python(std::ostream& out, PythonExportContext& ctx) override;
  virtual std::string export_vec_initial_partition_origin(std::ostream& out, PythonExportContext& ctx, const std::string& parent_name);


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

  double time_step;
  virtual void set_time_step(const double new_time_step_) {
    if (initialized) {
      throw RuntimeError("Value 'time_step' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    time_step = new_time_step_;
  }
  virtual double get_time_step() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return time_step;
  }

  double surface_grid_density;
  virtual void set_surface_grid_density(const double new_surface_grid_density_) {
    if (initialized) {
      throw RuntimeError("Value 'surface_grid_density' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    surface_grid_density = new_surface_grid_density_;
  }
  virtual double get_surface_grid_density() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return surface_grid_density;
  }

  double interaction_radius;
  virtual void set_interaction_radius(const double new_interaction_radius_) {
    if (initialized) {
      throw RuntimeError("Value 'interaction_radius' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    interaction_radius = new_interaction_radius_;
  }
  virtual double get_interaction_radius() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return interaction_radius;
  }

  double intermembrane_interaction_radius;
  virtual void set_intermembrane_interaction_radius(const double new_intermembrane_interaction_radius_) {
    if (initialized) {
      throw RuntimeError("Value 'intermembrane_interaction_radius' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    intermembrane_interaction_radius = new_intermembrane_interaction_radius_;
  }
  virtual double get_intermembrane_interaction_radius() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return intermembrane_interaction_radius;
  }

  double vacancy_search_distance;
  virtual void set_vacancy_search_distance(const double new_vacancy_search_distance_) {
    if (initialized) {
      throw RuntimeError("Value 'vacancy_search_distance' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    vacancy_search_distance = new_vacancy_search_distance_;
  }
  virtual double get_vacancy_search_distance() const {
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

  double partition_dimension;
  virtual void set_partition_dimension(const double new_partition_dimension_) {
    if (initialized) {
      throw RuntimeError("Value 'partition_dimension' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    partition_dimension = new_partition_dimension_;
  }
  virtual double get_partition_dimension() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return partition_dimension;
  }

  std::vector<double> initial_partition_origin;
  virtual void set_initial_partition_origin(const std::vector<double> new_initial_partition_origin_) {
    if (initialized) {
      throw RuntimeError("Value 'initial_partition_origin' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    initial_partition_origin = new_initial_partition_origin_;
  }
  virtual std::vector<double>& get_initial_partition_origin() {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return initial_partition_origin;
  }

  double subpartition_dimension;
  virtual void set_subpartition_dimension(const double new_subpartition_dimension_) {
    if (initialized) {
      throw RuntimeError("Value 'subpartition_dimension' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    subpartition_dimension = new_subpartition_dimension_;
  }
  virtual double get_subpartition_dimension() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return subpartition_dimension;
  }

  double total_iterations;
  virtual void set_total_iterations(const double new_total_iterations_) {
    if (initialized) {
      throw RuntimeError("Value 'total_iterations' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    total_iterations = new_total_iterations_;
  }
  virtual double get_total_iterations() const {
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

  int molecules_order_random_shuffle_periodicity;
  virtual void set_molecules_order_random_shuffle_periodicity(const int new_molecules_order_random_shuffle_periodicity_) {
    if (initialized) {
      throw RuntimeError("Value 'molecules_order_random_shuffle_periodicity' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    molecules_order_random_shuffle_periodicity = new_molecules_order_random_shuffle_periodicity_;
  }
  virtual int get_molecules_order_random_shuffle_periodicity() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return molecules_order_random_shuffle_periodicity;
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

  double initial_time;
  virtual void set_initial_time(const double new_initial_time_) {
    if (initialized) {
      throw RuntimeError("Value 'initial_time' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    initial_time = new_initial_time_;
  }
  virtual double get_initial_time() const {
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

  bool continue_after_sigalrm;
  virtual void set_continue_after_sigalrm(const bool new_continue_after_sigalrm_) {
    if (initialized) {
      throw RuntimeError("Value 'continue_after_sigalrm' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    continue_after_sigalrm = new_continue_after_sigalrm_;
  }
  virtual bool get_continue_after_sigalrm() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return continue_after_sigalrm;
  }

  // --- methods ---
}; // GenConfig

class Config;
py::class_<Config> define_pybinding_Config(py::module& m);
} // namespace API
} // namespace MCell

#endif // API_GEN_CONFIG_H
