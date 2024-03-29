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

#ifndef API_GEN_NOTIFICATIONS_H
#define API_GEN_NOTIFICATIONS_H

#include "api/api_common.h"
#include "api/base_data_class.h"

namespace MCell {
namespace API {

class Notifications;
class PythonExportContext;

#define NOTIFICATIONS_CTOR() \
    Notifications( \
        const int bng_verbosity_level_ = 0, \
        const bool rxn_and_species_report_ = false, \
        const int simulation_stats_every_n_iterations_ = 0, \
        const bool rxn_probability_changed_ = true, \
        const bool iteration_report_ = true, \
        const bool wall_overlap_report_ = false \
    ) { \
      class_name = "Notifications"; \
      bng_verbosity_level = bng_verbosity_level_; \
      rxn_and_species_report = rxn_and_species_report_; \
      simulation_stats_every_n_iterations = simulation_stats_every_n_iterations_; \
      rxn_probability_changed = rxn_probability_changed_; \
      iteration_report = iteration_report_; \
      wall_overlap_report = wall_overlap_report_; \
      postprocess_in_ctor(); \
      check_semantics(); \
    } \
    Notifications(DefaultCtorArgType) : \
      GenNotifications(DefaultCtorArgType()) { \
      set_all_attributes_as_default_or_unset(); \
      set_all_custom_attributes_to_default(); \
    }

class GenNotifications: public BaseDataClass {
public:
  GenNotifications() {
  }
  GenNotifications(DefaultCtorArgType) {
  }
  void postprocess_in_ctor() override {}
  void check_semantics() const override;
  void set_initialized() override;
  void set_all_attributes_as_default_or_unset() override;

  std::shared_ptr<Notifications> copy_notifications() const;
  std::shared_ptr<Notifications> deepcopy_notifications(py::dict = py::dict()) const;
  virtual bool __eq__(const Notifications& other) const;
  virtual bool eq_nonarray_attributes(const Notifications& other, const bool ignore_name = false) const;
  bool operator == (const Notifications& other) const { return __eq__(other);}
  bool operator != (const Notifications& other) const { return !__eq__(other);}
  std::string to_str(const bool all_details=false, const std::string ind="") const override;

  std::string export_to_python(std::ostream& out, PythonExportContext& ctx) override;


  // --- attributes ---
  int bng_verbosity_level;
  virtual void set_bng_verbosity_level(const int new_bng_verbosity_level_) {
    if (initialized) {
      throw RuntimeError("Value 'bng_verbosity_level' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    bng_verbosity_level = new_bng_verbosity_level_;
  }
  virtual int get_bng_verbosity_level() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return bng_verbosity_level;
  }

  bool rxn_and_species_report;
  virtual void set_rxn_and_species_report(const bool new_rxn_and_species_report_) {
    if (initialized) {
      throw RuntimeError("Value 'rxn_and_species_report' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    rxn_and_species_report = new_rxn_and_species_report_;
  }
  virtual bool get_rxn_and_species_report() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return rxn_and_species_report;
  }

  int simulation_stats_every_n_iterations;
  virtual void set_simulation_stats_every_n_iterations(const int new_simulation_stats_every_n_iterations_) {
    if (initialized) {
      throw RuntimeError("Value 'simulation_stats_every_n_iterations' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    simulation_stats_every_n_iterations = new_simulation_stats_every_n_iterations_;
  }
  virtual int get_simulation_stats_every_n_iterations() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return simulation_stats_every_n_iterations;
  }

  bool rxn_probability_changed;
  virtual void set_rxn_probability_changed(const bool new_rxn_probability_changed_) {
    if (initialized) {
      throw RuntimeError("Value 'rxn_probability_changed' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    rxn_probability_changed = new_rxn_probability_changed_;
  }
  virtual bool get_rxn_probability_changed() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return rxn_probability_changed;
  }

  bool iteration_report;
  virtual void set_iteration_report(const bool new_iteration_report_) {
    if (initialized) {
      throw RuntimeError("Value 'iteration_report' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    iteration_report = new_iteration_report_;
  }
  virtual bool get_iteration_report() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return iteration_report;
  }

  bool wall_overlap_report;
  virtual void set_wall_overlap_report(const bool new_wall_overlap_report_) {
    if (initialized) {
      throw RuntimeError("Value 'wall_overlap_report' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    wall_overlap_report = new_wall_overlap_report_;
  }
  virtual bool get_wall_overlap_report() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return wall_overlap_report;
  }

  // --- methods ---
}; // GenNotifications

class Notifications;
py::class_<Notifications> define_pybinding_Notifications(py::module& m);
} // namespace API
} // namespace MCell

#endif // API_GEN_NOTIFICATIONS_H
