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

#include "api/instantiation.h"

#include <algorithm>
#include "bng/bng.h"

#include "api/subsystem.h"
#include "api/geometry_object.h"
#include "api/surface_region.h"
#include "api/region.h"
#include "api/complex.h"

#include "generated/gen_geometry_utils.h"

using namespace std;

namespace MCell {
namespace API {

const char* const RELEASE_SITE_PREFIX = "rel_";


void Instantiation::dump() const {
  std::cout << to_str() << "\n";
}


std::shared_ptr<GeometryObject> Instantiation::find_volume_compartment(const std::string& name) {
  for (auto o: geometry_objects) {
    if (o->is_bngl_compartment) {
      if (o->name == name) {
        return o;
      }
    }
  }
  return std::shared_ptr<GeometryObject>(nullptr);
}


std::shared_ptr<GeometryObject> Instantiation::find_surface_compartment(const std::string& name) {
  for (auto o: geometry_objects) {
    if (o->is_bngl_compartment) {
      if (is_set(o->surface_compartment_name) && o->surface_compartment_name == name) {
        return o;
      }
    }
  }
  return std::shared_ptr<GeometryObject>(nullptr);
}


void Instantiation::load_bngl_seed_species(
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
  convert_bng_data_to_instantiation(bng_data, *subsystem, default_release_region);
}


void Instantiation::convert_bng_data_to_instantiation(
    const BNG::BNGData& bng_data,
    Subsystem& subsystem,
    std::shared_ptr<Region> default_release_region) {

  convert_compartments(bng_data);

  for (const BNG::SeedSpecies& bng_ss: bng_data.get_seed_species()) {
    convert_single_seed_species_to_release_site(bng_data, bng_ss, subsystem, default_release_region);
  }
}


void Instantiation::convert_compartments(const BNG::BNGData& bng_data) {

  if (bng_data.get_compartments().empty()) {
    return;
  }

  // create objects and assign 2d release regions
  for (const BNG::Compartment& bng_comp: bng_data.get_compartments()) {
    if (bng_comp.children_compartments.size() > 2) {
      throw RuntimeError("Automatic creation of compartments with more than one child compartment "
          "is not supported yet, error for '" + bng_comp.name + ".");
    }

    if (bng_comp.is_3d) {
      std::shared_ptr<GeometryObject> existing_obj = find_geometry_object(bng_comp.name);
      if (is_set(existing_obj)) {
        // object with this name already exists, we will use it instead
        existing_obj->is_bngl_compartment = true;
        continue;
      }

      float_t side = pow_f(bng_comp.get_volume(), 1.0/3.0);

      // create box for the given compartment
      shared_ptr<GeometryObject> box = geometry_utils::create_box(bng_comp.name, side);

      box->is_bngl_compartment = true;

      // set its 2d name
      if (bng_comp.has_parent()) {
        box->surface_compartment_name = bng_data.get_compartment(bng_comp.parent_compartment_id).name;
      }

      add_geometry_object(box);
    }
  }
}


void Instantiation::convert_single_seed_species_to_release_site(
    const BNG::BNGData& bng_data,
    const BNG::SeedSpecies& bng_ss,
    Subsystem& subsystem,
    std::shared_ptr<Region> default_release_region) {

  auto rel_site = make_shared<API::ReleaseSite>();

  // we need to create API representation for the cplx instance we got
  rel_site->complex =
      Subsystem::convert_cplx(bng_data, bng_ss.cplx);

  // releases from BNGL always used the default orientation
  rel_site->complex->orientation = Orientation::DEFAULT;

  if (bng_ss.cplx.has_compartment()) {
    bool surf_release = rel_site->complex->is_surf();
    const BNG::Compartment& c = bng_data.get_compartment(bng_ss.cplx.get_compartment_id());
    // check that dimensionality of compartment matches the released molecule
    if (surf_release && c.is_3d) {
      throw ValueError(S("Seed species specification for complex instance ") +
          bng_ss.cplx.to_str() + ": cannot release surface molecules " +
          "into a 3d compartment " + c.name + ".\n"
      );
    }
    else if (!surf_release && !c.is_3d) {
      throw ValueError(S("Seed species specification for complex instance ") +
          bng_ss.cplx.to_str() + ": cannot release volume molecules " +
          "onto a 2d compartment " + c.name + ".\n"
      );
    }

    rel_site->complex->compartment_name = c.name;
    rel_site->shape = Shape::COMPARTMENT;

    rel_site->name =
          "Release of " + rel_site->complex->to_bngl_str();
    if (is_set(rel_site->complex->compartment_name)) {
      rel_site->name += " at " + rel_site->complex->compartment_name;
    }
  }
  else {
    if (!is_set(default_release_region)) {
      throw ValueError(S("Seed species specification for complex instance ") +
          rel_site->complex->to_bngl_str() + " does not have a compartment and " +
          NAME_DEFAULT_RELEASE_REGION + " was set, don't know where to release.\n"
      );
    }

    rel_site->region = default_release_region;
    rel_site->shape = Shape::REGION_EXPR;
    rel_site->name =
          "Release of " + rel_site->complex->to_bngl_str() +
          " at " + rel_site->region->name;
  }


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


std::shared_ptr<Region> Instantiation::get_compartment_region(const std::string& name) {
  std::shared_ptr<GeometryObject> obj;

  // first try if it is a volume compartment
  obj = find_volume_compartment(name);
  if (is_set(obj)) {

    std::shared_ptr<Region> res = std::dynamic_pointer_cast<Region>(obj);

    for (auto child: obj->child_compartments) {
      res = res->__sub__(std::dynamic_pointer_cast<Region>(child));
    }
    return res;
  }

  // then a surface compartment
  obj = find_surface_compartment(name);
  if (is_set(obj)) {
    return std::dynamic_pointer_cast<Region>(obj);
  }

  // no compartment found
  return std::shared_ptr<Region>(nullptr);
}

}
// namespace API
} // namespace MCell
