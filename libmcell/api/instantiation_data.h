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

#ifndef API_INSTANTIATION_DATA_H
#define API_INSTANTIATION_DATA_H

#include "generated/gen_instantiation_data.h"
#include "api/common.h"
#include "api/api_utils.h"
#include "api/release_site.h"
#include "bng/bng_defines.h"

namespace BNG {
class BNGData;
class SeedSpecies;
class Compartment;
}

namespace MCell {
namespace API {

class Subsystem;

class InstantiationData: public GenInstantiationData {
public:

  typedef std::map<BNG::component_type_id_t, std::shared_ptr<Region>> CompartmentRegionMap;

  // from generated template
  void add_release_site(std::shared_ptr<ReleaseSite> s) override {
    append_to_vec(release_sites, s);
  }

  std::shared_ptr<ReleaseSite> find_release_site(const std::string& name) override {
    return vec_find_by_name(release_sites, name);
  }

  void add_geometry_object(std::shared_ptr<GeometryObject> o) override {
    append_to_vec(geometry_objects, o);
  }

  std::shared_ptr<GeometryObject> find_geometry_object(const std::string& name) override {
    return vec_find_by_name(geometry_objects, name);
  }

  void load_bngl_seed_species(
      const std::string& file_name,
      std::shared_ptr<Subsystem> subsystem,
      std::shared_ptr<Region> default_release_region = nullptr,
      const std::map<std::string, float_t>& parameter_overrides = std::map<std::string, float_t>()
  ) override;

  // added manually
  void dump() const;

protected:
  void convert_bng_data_to_instantiation_data(
      const BNG::BNGData& bng_data,
      Subsystem& subsystem,
      std::shared_ptr<Region> default_release_region);

private:
  void convert_compartments(
      const BNG::BNGData& bng_data,
      CompartmentRegionMap& compartment_region_map);

  void convert_single_seed_species_to_release_site(
      const BNG::BNGData& bng_data,
      const BNG::SeedSpecies& bng_ss,
      Subsystem& subsystem,
      std::shared_ptr<Region> default_release_region,
      const CompartmentRegionMap& compartment_region_map);
};

} // namespace API
} // namespace MCell

#endif // API_INSTANTIATION_DATA_H
