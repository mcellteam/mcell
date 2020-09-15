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

#include "api/instantiation_data.h"

#include "bng/bng.h"

#include "api/subsystem.h"
#include "api/complex_instance.h"

using namespace std;

namespace MCell {
namespace API {

const char* const RELEASE_SITE_PREFIX = "rel_";


void InstantiationData::dump() const {
  std::cout << to_str() << "\n";
}


void InstantiationData::load_bngl_seed_species(
    const std::string& file_name,
    std::shared_ptr<Subsystem> subsystem,
    std::shared_ptr<Region> default_release_region,
    const std::map<std::string, float_t>& parameter_overrides) {

  BNG::BNGData bng_data;

  int num_errors = BNG::parse_bngl_file(file_name, bng_data, parameter_overrides);
  if (num_errors != 0) {
    throw RuntimeError("Could not parse BNGL file " + file_name + ".");
  }

  // now convert everything we parsed into the API classes so that the user can
  // inspect or manipulate it if needed
  convert_bng_data_to_instantiation_data(bng_data, *subsystem, default_release_region);
}


void InstantiationData::convert_bng_data_to_instantiation_data(
    const BNG::BNGData& bng_data,
    Subsystem& subsystem,
    std::shared_ptr<Region> default_release_region) {

  for (const BNG::SeedSpecies& bng_ss: bng_data.get_seed_species()) {
    convert_single_seed_species_to_release_site(bng_data, bng_ss, subsystem, default_release_region);
  }
}

void InstantiationData::convert_single_seed_species_to_release_site(
    const BNG::BNGData& bng_data,
    const BNG::SeedSpecies& bng_ss,
    Subsystem& subsystem,
    std::shared_ptr<Region> default_release_region) {

  if (!is_set(default_release_region)) {
    throw ValueError(S("Parameter ") + NAME_DEFAULT_RELEASE_REGION +
        " must be currently always set because compartments are not supported yet.");
  }

  auto rel_site = make_shared<API::ReleaseSite>();

  // we need to create API representation for the cplx instance we got
  rel_site->complex_instance =
      subsystem.convert_cplx_instance(bng_data, bng_ss.cplx_instance);

  rel_site->name =
      "Release of " + rel_site->complex_instance->to_bngl_str() +
      " at " + default_release_region->name;

  rel_site->region = default_release_region;

  uint truncated_count = BNG::floor_f(bng_ss.count);
  if (bng_ss.count != truncated_count) {
    cout << "Warning: Release count of complex instance created from loaded BNGL file '" +
        rel_site->name + "' was truncated from " << bng_ss.count << " to " << truncated_count << ".\n";
  }
  rel_site->number_to_release = truncated_count;

  rel_site->check_semantics(); // only for internal checks
  rel_site->postprocess_in_ctor(); // sets shape

  release_sites.push_back(rel_site);
}

}
// namespace API
} // namespace MCell
