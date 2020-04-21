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
        const std::string& shape_, \
        std::shared_ptr<Species> molecule_, \
        const Vec3& location_ = VEC3_UNSET, \
        const float_t site_diameter_ = FLT_UNSET, \
        const float_t site_radius_ = FLT_UNSET, \
        const float_t release_probability_ = FLT_UNSET \
    ) { \
      class_name = "ReleaseSite"; \
      name = name_; \
      shape = shape_; \
      molecule = molecule_; \
      location = location_; \
      site_diameter = site_diameter_; \
      site_radius = site_radius_; \
      release_probability = release_probability_; \
    }

class GenReleaseSite: public BaseDataClass {
public:
  SemRes check_semantics(std::ostream& out) const override;
  std::string to_str() const override;

  // --- attributes ---
  std::string shape;
  virtual void set_shape(const std::string& new_shape_) {
    shape = new_shape_;
  }
  virtual const std::string& get_shape() const {
    return shape;
  }

  std::shared_ptr<Species> molecule;
  virtual void set_molecule(std::shared_ptr<Species> new_molecule_) {
    molecule = new_molecule_;
  }
  virtual std::shared_ptr<Species> get_molecule() const {
    return molecule;
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
