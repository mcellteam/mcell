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

#ifndef API_GEN_NOTIFICATIONS_H
#define API_GEN_NOTIFICATIONS_H

#include "api/common.h"
#include "api/base_data_class.h"

namespace MCell {
namespace API {

class Notifications;
class PythonExportContext;

#define NOTIFICATIONS_CTOR() \
    Notifications( \
        const int bng_verbosity_level_ = 0, \
        const bool rxn_and_species_report_ = true, \
        const int simulation_stats_every_n_iterations_ = 0 \
    ) { \
      class_name = "Notifications"; \
      bng_verbosity_level = bng_verbosity_level_; \
      rxn_and_species_report = rxn_and_species_report_; \
      simulation_stats_every_n_iterations = simulation_stats_every_n_iterations_; \
      postprocess_in_ctor();\
      check_semantics();\
    }

class GenNotifications: public BaseDataClass {
public:
  void postprocess_in_ctor() override {}
  void check_semantics() const override;
  void set_initialized() override;
  void set_all_attributes_as_default_or_unset() override;

  virtual bool __eq__(const Notifications& other) const;
  virtual bool eq_nonarray_attributes(const Notifications& other, const bool ignore_name = false) const;
  bool operator == (const Notifications& other) const { return __eq__(other);}
  bool operator != (const Notifications& other) const { return !__eq__(other);}
  std::string to_str(const std::string ind="") const override;

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

  // --- methods ---
}; // GenNotifications

class Notifications;
py::class_<Notifications> define_pybinding_Notifications(py::module& m);
} // namespace API
} // namespace MCell

#endif // API_GEN_NOTIFICATIONS_H
