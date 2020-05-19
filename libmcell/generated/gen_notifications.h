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

#include "../api/common.h"
#include "../api/base_data_class.h"

namespace MCell {
namespace API {

#define NOTIFICATIONS_CTOR() \
    Notifications( \
        const bool probability_report_ = true, \
        const Notification diffusion_constant_report_ = Notification::Brief, \
        const bool final_summary_ = true, \
        const bool iteration_report_ = true, \
        const bool varying_probability_report_ = true, \
        const bool progress_report_ = true, \
        const bool release_event_report_ = true, \
        const bool molecule_collision_report_ = true \
    ) { \
      class_name = "Notifications"; \
      probability_report = probability_report_; \
      diffusion_constant_report = diffusion_constant_report_; \
      final_summary = final_summary_; \
      iteration_report = iteration_report_; \
      varying_probability_report = varying_probability_report_; \
      progress_report = progress_report_; \
      release_event_report = release_event_report_; \
      molecule_collision_report = molecule_collision_report_; \
      postprocess_in_ctor();\
      check_semantics();\
    }

class GenNotifications: public BaseDataClass {
public:
  void postprocess_in_ctor() override {}
  void check_semantics() const override;
  std::string to_str(const std::string ind="") const override;

  void set_initialized() override;

  bool __eq__(const GenNotifications& other) const;
  // --- attributes ---
  bool probability_report;
  virtual void set_probability_report(const bool new_probability_report_) {
    if (initialized) {
      throw RuntimeError("Value 'probability_report' of object with name " + name + " (class " + class_name + ")"
                         "cannot be set after model was initialized.");
    }
    probability_report = new_probability_report_;
  }
  virtual bool get_probability_report() const {
    return probability_report;
  }

  Notification diffusion_constant_report;
  virtual void set_diffusion_constant_report(const Notification new_diffusion_constant_report_) {
    if (initialized) {
      throw RuntimeError("Value 'diffusion_constant_report' of object with name " + name + " (class " + class_name + ")"
                         "cannot be set after model was initialized.");
    }
    diffusion_constant_report = new_diffusion_constant_report_;
  }
  virtual Notification get_diffusion_constant_report() const {
    return diffusion_constant_report;
  }

  bool final_summary;
  virtual void set_final_summary(const bool new_final_summary_) {
    if (initialized) {
      throw RuntimeError("Value 'final_summary' of object with name " + name + " (class " + class_name + ")"
                         "cannot be set after model was initialized.");
    }
    final_summary = new_final_summary_;
  }
  virtual bool get_final_summary() const {
    return final_summary;
  }

  bool iteration_report;
  virtual void set_iteration_report(const bool new_iteration_report_) {
    if (initialized) {
      throw RuntimeError("Value 'iteration_report' of object with name " + name + " (class " + class_name + ")"
                         "cannot be set after model was initialized.");
    }
    iteration_report = new_iteration_report_;
  }
  virtual bool get_iteration_report() const {
    return iteration_report;
  }

  bool varying_probability_report;
  virtual void set_varying_probability_report(const bool new_varying_probability_report_) {
    if (initialized) {
      throw RuntimeError("Value 'varying_probability_report' of object with name " + name + " (class " + class_name + ")"
                         "cannot be set after model was initialized.");
    }
    varying_probability_report = new_varying_probability_report_;
  }
  virtual bool get_varying_probability_report() const {
    return varying_probability_report;
  }

  bool progress_report;
  virtual void set_progress_report(const bool new_progress_report_) {
    if (initialized) {
      throw RuntimeError("Value 'progress_report' of object with name " + name + " (class " + class_name + ")"
                         "cannot be set after model was initialized.");
    }
    progress_report = new_progress_report_;
  }
  virtual bool get_progress_report() const {
    return progress_report;
  }

  bool release_event_report;
  virtual void set_release_event_report(const bool new_release_event_report_) {
    if (initialized) {
      throw RuntimeError("Value 'release_event_report' of object with name " + name + " (class " + class_name + ")"
                         "cannot be set after model was initialized.");
    }
    release_event_report = new_release_event_report_;
  }
  virtual bool get_release_event_report() const {
    return release_event_report;
  }

  bool molecule_collision_report;
  virtual void set_molecule_collision_report(const bool new_molecule_collision_report_) {
    if (initialized) {
      throw RuntimeError("Value 'molecule_collision_report' of object with name " + name + " (class " + class_name + ")"
                         "cannot be set after model was initialized.");
    }
    molecule_collision_report = new_molecule_collision_report_;
  }
  virtual bool get_molecule_collision_report() const {
    return molecule_collision_report;
  }

  // --- methods ---
}; // GenNotifications

class Notifications;
py::class_<Notifications> define_pybinding_Notifications(py::module& m);
} // namespace API
} // namespace MCell

#endif // API_GEN_NOTIFICATIONS_H
