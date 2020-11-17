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

#ifndef API_GEN_INSTANTIATION_DATA_H
#define API_GEN_INSTANTIATION_DATA_H

#include "api/common.h"

namespace MCell {
namespace API {

class InstantiationData;
class GeometryObject;
class Region;
class ReleaseSite;
class Subsystem;

class GenInstantiationData {
public:
  virtual ~GenInstantiationData() {}
  virtual bool __eq__(const InstantiationData& other) const;
  bool operator == (const InstantiationData& other) const { return __eq__(other);}
  bool operator != (const InstantiationData& other) const { return !__eq__(other);}
  std::string to_str(const std::string ind="") const ;

  // --- attributes ---
  std::vector<std::shared_ptr<ReleaseSite>> release_sites;
  virtual void set_release_sites(const std::vector<std::shared_ptr<ReleaseSite>> new_release_sites_) {
    release_sites = new_release_sites_;
  }
  virtual std::vector<std::shared_ptr<ReleaseSite>> get_release_sites() const {
    return release_sites;
  }

  std::vector<std::shared_ptr<GeometryObject>> geometry_objects;
  virtual void set_geometry_objects(const std::vector<std::shared_ptr<GeometryObject>> new_geometry_objects_) {
    geometry_objects = new_geometry_objects_;
  }
  virtual std::vector<std::shared_ptr<GeometryObject>> get_geometry_objects() const {
    return geometry_objects;
  }

  // --- methods ---
  virtual void add_release_site(std::shared_ptr<ReleaseSite> s) = 0;
  virtual std::shared_ptr<ReleaseSite> find_release_site(const std::string& name) = 0;
  virtual void add_geometry_object(std::shared_ptr<GeometryObject> o) = 0;
  virtual std::shared_ptr<GeometryObject> find_geometry_object(const std::string& name) = 0;
  virtual std::shared_ptr<GeometryObject> find_volume_compartment(const std::string& name) = 0;
  virtual std::shared_ptr<GeometryObject> find_surface_compartment(const std::string& name) = 0;
  virtual void load_bngl_seed_species(const std::string& file_name, std::shared_ptr<Subsystem> subsystem, std::shared_ptr<Region> default_release_region = nullptr, const std::map<std::string, float_t>& parameter_overrides = std::map<std::string, float_t>()) = 0;
}; // GenInstantiationData

class InstantiationData;
py::class_<InstantiationData> define_pybinding_InstantiationData(py::module& m);
} // namespace API
} // namespace MCell

#endif // API_GEN_INSTANTIATION_DATA_H
