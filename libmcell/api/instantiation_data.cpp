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

#include "generated/gen_geometry_utils.h"

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

  CompartmentRegionMap compartment_region_map;
  convert_compartments(bng_data, compartment_region_map);

  for (const BNG::SeedSpecies& bng_ss: bng_data.get_seed_species()) {
    convert_single_seed_species_to_release_site(
        bng_data, bng_ss, subsystem, default_release_region, compartment_region_map);
  }
}


void InstantiationData::convert_compartments(
    const BNG::BNGData& bng_data,
    CompartmentRegionMap& compartment_region_map) {

  // check that each parent has a single child and that
  // the lowest level child is 3d, its parent must be 2d,
  // its 3d, etc.
  // we need to find the one that is no parent to any compartment
  // map<parent, child>
  map<string, string> hierarchy;
  for (const BNG::Compartment& c: bng_data.get_compartments()) {
    if (c.parent_name != "") {
      // check that parent has up to 1 child
      if (hierarchy.count(c.parent_name) != 0) {
        throw RuntimeError(
            "Single compartment may have max one child, error for parent compartment " + c.parent_name + ".");
      }
      // check that 2d/3d changes
      const BNG::Compartment& parent =
          bng_data.get_compartment(bng_data.find_compartment_id(c.parent_name));
      if (c.is_3d == parent.is_3d) {
        // NOTE: this check should be rather in BNG lib
        throw RuntimeError(
            "Parent compartment " + parent.name +
            " must be of different dimensionality than its child " + c.name + ".");
      }

      hierarchy[c.parent_name] = c.name;
    }
  }


  // create objects
  for (const BNG::Compartment& bng_comp: bng_data.get_compartments()) {
    if (bng_comp.is_3d) {
      float_t side = pow_f(bng_comp.volume, 1.0/3.0);
      std::shared_ptr<GeometryObject> box = geometry_utils::create_box(bng_comp.name, side);
      geometry_objects.push_back(box);

      // 3d compartment corresponds to the created box
      compartment_region_map[bng_comp.name] = box;

      // the same goes for its 2d parent
      if (bng_comp.parent_name != "") {
        compartment_region_map[bng_comp.parent_name] = box;
      }
    }
  }
}


void InstantiationData::convert_single_seed_species_to_release_site(
    const BNG::BNGData& bng_data,
    const BNG::SeedSpecies& bng_ss,
    Subsystem& subsystem,
    std::shared_ptr<Region> default_release_region,
    const CompartmentRegionMap& compartment_region_map) {

  auto rel_site = make_shared<API::ReleaseSite>();

  // we need to create API representation for the cplx instance we got
  rel_site->complex_instance =
      subsystem.convert_cplx_instance(bng_data, bng_ss.cplx_instance);

  bool surf_release = rel_site->complex_instance->is_surf();
  if (surf_release) {
    // the default orientation of released molecules is 'up'
    rel_site->orientation = Orientation::UP;
  }

  if (bng_ss.compartment_id != BNG::COMPARTMENT_ID_INVALID) {
    const BNG::Compartment& c = bng_data.get_compartment(bng_ss.compartment_id);
    // check that dimensionality of compartment matches the released molecule
    if (surf_release && c.is_3d) {
      throw ValueError(S("Seed species specification for complex instance ") +
          rel_site->complex_instance->name + ": cannot release surface molecules " +
          "into a 3d compartment " + c.name + ".\n"
      );
    }
    else if (!surf_release && !c.is_3d) {
      throw ValueError(S("Seed species specification for complex instance ") +
          rel_site->complex_instance->name + ": cannot release volume molecules " +
          "onto a 2d compartment " + c.name + ".\n"
      );
    }

    auto it = compartment_region_map.find(c.name);
    assert(it != compartment_region_map.end());
    rel_site->region = it->second;
  }
  else {
    if (!is_set(default_release_region)) {
      throw ValueError(S("Seed species specification for complex instance ") +
          rel_site->complex_instance->name + " does not have a compartment and neither " +
          NAME_DEFAULT_RELEASE_REGION + " was set, don't know where to release.\n"
      );
    }

    rel_site->region = default_release_region;
  }

  rel_site->name =
      "Release of " + rel_site->complex_instance->to_bngl_str() +
      " at " + rel_site->region->name;

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
