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

#ifndef API_GEN_WARNINGS_H
#define API_GEN_WARNINGS_H

#include "../api/common.h"
#include "../api/base_data_class.h"

namespace MCell {
namespace API {

#define WARNINGS_CTOR() \
    Warnings( \
        const WarningLevel molecule_collision_report_ = WarningLevel::WARNING, \
        const WarningLevel degenerate_polygons_ = WarningLevel::WARNING, \
        const WarningLevel negative_diffusion_constant_ = WarningLevel::WARNING, \
        const WarningLevel missing_surface_orientation_ = WarningLevel::ERROR, \
        const WarningLevel negative_reaction_rate_ = WarningLevel::WARNING, \
        const WarningLevel useless_volume_orientation_ = WarningLevel::WARNING, \
        const WarningLevel high_reaction_probability_ = WarningLevel::IGNORE, \
        const WarningLevel lifetime_too_short_ = WarningLevel::WARNING, \
        const float_t lifetime_threshold_ = 50, \
        const WarningLevel missed_reactions_ = WarningLevel::WARNING, \
        const float_t missed_reactions_threshold_ = 0.00100000004749745 \
    ) { \
      class_name = "Warnings"; \
      molecule_collision_report = molecule_collision_report_; \
      degenerate_polygons = degenerate_polygons_; \
      negative_diffusion_constant = negative_diffusion_constant_; \
      missing_surface_orientation = missing_surface_orientation_; \
      negative_reaction_rate = negative_reaction_rate_; \
      useless_volume_orientation = useless_volume_orientation_; \
      high_reaction_probability = high_reaction_probability_; \
      lifetime_too_short = lifetime_too_short_; \
      lifetime_threshold = lifetime_threshold_; \
      missed_reactions = missed_reactions_; \
      missed_reactions_threshold = missed_reactions_threshold_; \
      postprocess_in_ctor();\
      check_semantics();\
    }

class GenWarnings: public BaseDataClass {
public:
  void postprocess_in_ctor() override {}
  void check_semantics() const override;
  void set_initialized() override;
  bool __eq__(const GenWarnings& other) const;
  void set_all_attributes_as_default_or_unset() override;

  std::string to_str(const std::string ind="") const override;

  // --- attributes ---
  WarningLevel molecule_collision_report;
  virtual void set_molecule_collision_report(const WarningLevel new_molecule_collision_report_) {
    if (initialized) {
      throw RuntimeError("Value 'molecule_collision_report' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    molecule_collision_report = new_molecule_collision_report_;
  }
  virtual WarningLevel get_molecule_collision_report() const {
    return molecule_collision_report;
  }

  WarningLevel degenerate_polygons;
  virtual void set_degenerate_polygons(const WarningLevel new_degenerate_polygons_) {
    if (initialized) {
      throw RuntimeError("Value 'degenerate_polygons' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    degenerate_polygons = new_degenerate_polygons_;
  }
  virtual WarningLevel get_degenerate_polygons() const {
    return degenerate_polygons;
  }

  WarningLevel negative_diffusion_constant;
  virtual void set_negative_diffusion_constant(const WarningLevel new_negative_diffusion_constant_) {
    if (initialized) {
      throw RuntimeError("Value 'negative_diffusion_constant' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    negative_diffusion_constant = new_negative_diffusion_constant_;
  }
  virtual WarningLevel get_negative_diffusion_constant() const {
    return negative_diffusion_constant;
  }

  WarningLevel missing_surface_orientation;
  virtual void set_missing_surface_orientation(const WarningLevel new_missing_surface_orientation_) {
    if (initialized) {
      throw RuntimeError("Value 'missing_surface_orientation' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    missing_surface_orientation = new_missing_surface_orientation_;
  }
  virtual WarningLevel get_missing_surface_orientation() const {
    return missing_surface_orientation;
  }

  WarningLevel negative_reaction_rate;
  virtual void set_negative_reaction_rate(const WarningLevel new_negative_reaction_rate_) {
    if (initialized) {
      throw RuntimeError("Value 'negative_reaction_rate' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    negative_reaction_rate = new_negative_reaction_rate_;
  }
  virtual WarningLevel get_negative_reaction_rate() const {
    return negative_reaction_rate;
  }

  WarningLevel useless_volume_orientation;
  virtual void set_useless_volume_orientation(const WarningLevel new_useless_volume_orientation_) {
    if (initialized) {
      throw RuntimeError("Value 'useless_volume_orientation' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    useless_volume_orientation = new_useless_volume_orientation_;
  }
  virtual WarningLevel get_useless_volume_orientation() const {
    return useless_volume_orientation;
  }

  WarningLevel high_reaction_probability;
  virtual void set_high_reaction_probability(const WarningLevel new_high_reaction_probability_) {
    if (initialized) {
      throw RuntimeError("Value 'high_reaction_probability' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    high_reaction_probability = new_high_reaction_probability_;
  }
  virtual WarningLevel get_high_reaction_probability() const {
    return high_reaction_probability;
  }

  WarningLevel lifetime_too_short;
  virtual void set_lifetime_too_short(const WarningLevel new_lifetime_too_short_) {
    if (initialized) {
      throw RuntimeError("Value 'lifetime_too_short' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    lifetime_too_short = new_lifetime_too_short_;
  }
  virtual WarningLevel get_lifetime_too_short() const {
    return lifetime_too_short;
  }

  float_t lifetime_threshold;
  virtual void set_lifetime_threshold(const float_t new_lifetime_threshold_) {
    if (initialized) {
      throw RuntimeError("Value 'lifetime_threshold' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    lifetime_threshold = new_lifetime_threshold_;
  }
  virtual float_t get_lifetime_threshold() const {
    return lifetime_threshold;
  }

  WarningLevel missed_reactions;
  virtual void set_missed_reactions(const WarningLevel new_missed_reactions_) {
    if (initialized) {
      throw RuntimeError("Value 'missed_reactions' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    missed_reactions = new_missed_reactions_;
  }
  virtual WarningLevel get_missed_reactions() const {
    return missed_reactions;
  }

  float_t missed_reactions_threshold;
  virtual void set_missed_reactions_threshold(const float_t new_missed_reactions_threshold_) {
    if (initialized) {
      throw RuntimeError("Value 'missed_reactions_threshold' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    missed_reactions_threshold = new_missed_reactions_threshold_;
  }
  virtual float_t get_missed_reactions_threshold() const {
    return missed_reactions_threshold;
  }

  // --- methods ---
}; // GenWarnings

class Warnings;
py::class_<Warnings> define_pybinding_Warnings(py::module& m);
} // namespace API
} // namespace MCell

#endif // API_GEN_WARNINGS_H
