/******************************************************************************
 *
 * Copyright (C) 2020 by
 * The Salk Institute for Biological Studies
 *
 * Use of this source code is governed by an MIT-style
 * license that can be found in the LICENSE file or at
 * https://opensource.org/licenses/MIT.
 *
******************************************************************************/

#ifndef API_INSTANTIATION_H
#define API_INSTANTIATION_H

#include "generated/gen_instantiation.h"
#include "api/api_common.h"
#include "api/api_utils.h"
#include "api/release_site.h"

namespace BNG {
class BNGData;
class SeedSpecies;
class Compartment;
}

namespace MCell {
namespace API {

class Instantiation: public GenInstantiation {
public:
  INSTANTIATION_CTOR()

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

  void load_bngl_compartments_and_seed_species(
      const std::string& file_name,
      std::shared_ptr<Region> default_release_region = nullptr,
      const std::map<std::string, double>& parameter_overrides = std::map<std::string, double>()
  ) override;

  std::shared_ptr<GeometryObject> find_volume_compartment_object(const std::string& name) override;
  std::shared_ptr<GeometryObject> find_surface_compartment_object(const std::string& name) override;

  // added manually
  void dump() const;

  // returns empty shared ptr if compartment was not found
  std::shared_ptr<Region> get_compartment_region(const std::string& name);

protected:
  void convert_bng_data_to_instantiation(
      const BNG::BNGData& bng_data,
      std::shared_ptr<Region> default_release_region);

private:
  void convert_compartments(
      const BNG::BNGData& bng_data);

  void convert_single_seed_species_to_release_site(
      const BNG::BNGData& bng_data,
      const BNG::SeedSpecies& bng_ss,
      std::shared_ptr<Region> default_release_region);
};

} // namespace API
} // namespace MCell

#endif // API_INSTANTIATION_H
