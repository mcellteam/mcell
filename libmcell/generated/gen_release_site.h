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

#ifndef API_GEN_RELEASE_SITE_H
#define API_GEN_RELEASE_SITE_H

#include "../api/common.h"
#include "../api/base_data_class.h"

namespace MCell {
namespace API {

class Complex;
class MoleculeReleaseInfo;
class Region;
class ReleasePattern;

#define RELEASE_SITE_CTOR() \
    ReleaseSite( \
        const std::string& name_, \
        std::shared_ptr<Complex> complex_ = nullptr, \
        const Orientation orientation_ = Orientation::NONE, \
        const std::vector<std::shared_ptr<MoleculeReleaseInfo>> molecule_list_ = std::vector<std::shared_ptr<MoleculeReleaseInfo>>(), \
        const float_t release_time_ = 0, \
        std::shared_ptr<ReleasePattern> release_pattern_ = nullptr, \
        const Shape shape_ = Shape::UNSET, \
        std::shared_ptr<Region> region_ = nullptr, \
        const std::string& compartment_name_ = STR_UNSET, \
        const Vec3& location_ = VEC3_UNSET, \
        const float_t site_diameter_ = 0, \
        const float_t site_radius_ = FLT_UNSET, \
        const float_t number_to_release_ = FLT_UNSET, \
        const float_t density_ = FLT_UNSET, \
        const float_t concentration_ = FLT_UNSET, \
        const float_t release_probability_ = FLT_UNSET \
    ) { \
      class_name = "ReleaseSite"; \
      name = name_; \
      complex = complex_; \
      orientation = orientation_; \
      molecule_list = molecule_list_; \
      release_time = release_time_; \
      release_pattern = release_pattern_; \
      shape = shape_; \
      region = region_; \
      compartment_name = compartment_name_; \
      location = location_; \
      site_diameter = site_diameter_; \
      site_radius = site_radius_; \
      number_to_release = number_to_release_; \
      density = density_; \
      concentration = concentration_; \
      release_probability = release_probability_; \
      postprocess_in_ctor();\
      check_semantics();\
    }

class GenReleaseSite: public BaseDataClass {
public:
  void postprocess_in_ctor() override {}
  void check_semantics() const override;
  void set_initialized() override;
  bool __eq__(const GenReleaseSite& other) const;
  void set_all_attributes_as_default_or_unset() override;

  std::string to_str(const std::string ind="") const override;

  // --- attributes ---
  std::shared_ptr<Complex> complex;
  virtual void set_complex(std::shared_ptr<Complex> new_complex_) {
    if (initialized) {
      throw RuntimeError("Value 'complex' of object with name " + name + " (class " + class_name + ")"
                         "cannot be set after model was initialized.");
    }
    complex = new_complex_;
  }
  virtual std::shared_ptr<Complex> get_complex() const {
    return complex;
  }

  Orientation orientation;
  virtual void set_orientation(const Orientation new_orientation_) {
    if (initialized) {
      throw RuntimeError("Value 'orientation' of object with name " + name + " (class " + class_name + ")"
                         "cannot be set after model was initialized.");
    }
    orientation = new_orientation_;
  }
  virtual Orientation get_orientation() const {
    return orientation;
  }

  std::vector<std::shared_ptr<MoleculeReleaseInfo>> molecule_list;
  virtual void set_molecule_list(const std::vector<std::shared_ptr<MoleculeReleaseInfo>> new_molecule_list_) {
    if (initialized) {
      throw RuntimeError("Value 'molecule_list' of object with name " + name + " (class " + class_name + ")"
                         "cannot be set after model was initialized.");
    }
    molecule_list = new_molecule_list_;
  }
  virtual std::vector<std::shared_ptr<MoleculeReleaseInfo>> get_molecule_list() const {
    return molecule_list;
  }

  float_t release_time;
  virtual void set_release_time(const float_t new_release_time_) {
    if (initialized) {
      throw RuntimeError("Value 'release_time' of object with name " + name + " (class " + class_name + ")"
                         "cannot be set after model was initialized.");
    }
    release_time = new_release_time_;
  }
  virtual float_t get_release_time() const {
    return release_time;
  }

  std::shared_ptr<ReleasePattern> release_pattern;
  virtual void set_release_pattern(std::shared_ptr<ReleasePattern> new_release_pattern_) {
    if (initialized) {
      throw RuntimeError("Value 'release_pattern' of object with name " + name + " (class " + class_name + ")"
                         "cannot be set after model was initialized.");
    }
    release_pattern = new_release_pattern_;
  }
  virtual std::shared_ptr<ReleasePattern> get_release_pattern() const {
    return release_pattern;
  }

  Shape shape;
  virtual void set_shape(const Shape new_shape_) {
    if (initialized) {
      throw RuntimeError("Value 'shape' of object with name " + name + " (class " + class_name + ")"
                         "cannot be set after model was initialized.");
    }
    shape = new_shape_;
  }
  virtual Shape get_shape() const {
    return shape;
  }

  std::shared_ptr<Region> region;
  virtual void set_region(std::shared_ptr<Region> new_region_) {
    if (initialized) {
      throw RuntimeError("Value 'region' of object with name " + name + " (class " + class_name + ")"
                         "cannot be set after model was initialized.");
    }
    region = new_region_;
  }
  virtual std::shared_ptr<Region> get_region() const {
    return region;
  }

  std::string compartment_name;
  virtual void set_compartment_name(const std::string& new_compartment_name_) {
    if (initialized) {
      throw RuntimeError("Value 'compartment_name' of object with name " + name + " (class " + class_name + ")"
                         "cannot be set after model was initialized.");
    }
    compartment_name = new_compartment_name_;
  }
  virtual const std::string& get_compartment_name() const {
    return compartment_name;
  }

  Vec3 location;
  virtual void set_location(const Vec3& new_location_) {
    if (initialized) {
      throw RuntimeError("Value 'location' of object with name " + name + " (class " + class_name + ")"
                         "cannot be set after model was initialized.");
    }
    location = new_location_;
  }
  virtual const Vec3& get_location() const {
    return location;
  }

  float_t site_diameter;
  virtual void set_site_diameter(const float_t new_site_diameter_) {
    if (initialized) {
      throw RuntimeError("Value 'site_diameter' of object with name " + name + " (class " + class_name + ")"
                         "cannot be set after model was initialized.");
    }
    site_diameter = new_site_diameter_;
  }
  virtual float_t get_site_diameter() const {
    return site_diameter;
  }

  float_t site_radius;
  virtual void set_site_radius(const float_t new_site_radius_) {
    if (initialized) {
      throw RuntimeError("Value 'site_radius' of object with name " + name + " (class " + class_name + ")"
                         "cannot be set after model was initialized.");
    }
    site_radius = new_site_radius_;
  }
  virtual float_t get_site_radius() const {
    return site_radius;
  }

  float_t number_to_release;
  virtual void set_number_to_release(const float_t new_number_to_release_) {
    if (initialized) {
      throw RuntimeError("Value 'number_to_release' of object with name " + name + " (class " + class_name + ")"
                         "cannot be set after model was initialized.");
    }
    number_to_release = new_number_to_release_;
  }
  virtual float_t get_number_to_release() const {
    return number_to_release;
  }

  float_t density;
  virtual void set_density(const float_t new_density_) {
    if (initialized) {
      throw RuntimeError("Value 'density' of object with name " + name + " (class " + class_name + ")"
                         "cannot be set after model was initialized.");
    }
    density = new_density_;
  }
  virtual float_t get_density() const {
    return density;
  }

  float_t concentration;
  virtual void set_concentration(const float_t new_concentration_) {
    if (initialized) {
      throw RuntimeError("Value 'concentration' of object with name " + name + " (class " + class_name + ")"
                         "cannot be set after model was initialized.");
    }
    concentration = new_concentration_;
  }
  virtual float_t get_concentration() const {
    return concentration;
  }

  float_t release_probability;
  virtual void set_release_probability(const float_t new_release_probability_) {
    if (initialized) {
      throw RuntimeError("Value 'release_probability' of object with name " + name + " (class " + class_name + ")"
                         "cannot be set after model was initialized.");
    }
    release_probability = new_release_probability_;
  }
  virtual float_t get_release_probability() const {
    return release_probability;
  }

  // --- methods ---
}; // GenReleaseSite

class ReleaseSite;
py::class_<ReleaseSite> define_pybinding_ReleaseSite(py::module& m);
} // namespace API
} // namespace MCell

#endif // API_GEN_RELEASE_SITE_H
