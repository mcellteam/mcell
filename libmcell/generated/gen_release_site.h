/******************************************************************************
 *
 * Copyright (C) 2020 by
 * The Salk Institute for Biological Studies and
 * Pittsburgh Supercomputing Center, Carnegie Mellon University
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

namespace MCell {
namespace API {

class Species;

#define RELEASE_SITE_CTOR() \
    ReleaseSite( \
        const std::string& name_, \
        std::shared_ptr<Species> species_, \
        const Shape shape_ = Shape::Unset, \
        const Vec3& location_ = VEC3_UNSET, \
        const float_t site_diameter_ = FLT_UNSET, \
        const float_t site_radius_ = FLT_UNSET, \
        const int number_to_release_ = INT_UNSET, \
        const float_t release_probability_ = FLT_UNSET \
    ) { \
      class_name = "ReleaseSite"; \
      name = name_; \
      species = species_; \
      shape = shape_; \
      location = location_; \
      site_diameter = site_diameter_; \
      site_radius = site_radius_; \
      number_to_release = number_to_release_; \
      release_probability = release_probability_; \
      postprocess_in_ctor();\
      check_semantics();\
    }

class GenReleaseSite: public BaseDataClass {
public:
  void postprocess_in_ctor() override {}
  void check_semantics() const override;
  std::string to_str(const std::string ind="") const override;

  // --- attributes ---
  std::shared_ptr<Species> species;
  virtual void set_species(std::shared_ptr<Species> new_species_) {
    species = new_species_;
  }
  virtual std::shared_ptr<Species> get_species() const {
    return species;
  }

  Shape shape;
  virtual void set_shape(const Shape new_shape_) {
    shape = new_shape_;
  }
  virtual Shape get_shape() const {
    return shape;
  }

  Vec3 location;
  virtual void set_location(const Vec3& new_location_) {
    location = new_location_;
  }
  virtual const Vec3& get_location() const {
    return location;
  }

  float_t site_diameter;
  virtual void set_site_diameter(const float_t new_site_diameter_) {
    site_diameter = new_site_diameter_;
  }
  virtual float_t get_site_diameter() const {
    return site_diameter;
  }

  float_t site_radius;
  virtual void set_site_radius(const float_t new_site_radius_) {
    site_radius = new_site_radius_;
  }
  virtual float_t get_site_radius() const {
    return site_radius;
  }

  int number_to_release;
  virtual void set_number_to_release(const int new_number_to_release_) {
    number_to_release = new_number_to_release_;
  }
  virtual int get_number_to_release() const {
    return number_to_release;
  }

  float_t release_probability;
  virtual void set_release_probability(const float_t new_release_probability_) {
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
